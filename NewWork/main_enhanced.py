from itertools import permutations
import time
from collections import defaultdict
from pathlib import Path


class SBoxImplicitCalculator:
    """S盒隐式方程计算器 - 高斯消元法"""

    def __init__(self, n_bits=3):
        self.n_bits = n_bits
        self.n_inputs = 1 << n_bits

    def generate_all_monomials(self, max_degree=2):
        """生成所有可能的单项式（次数≤max_degree）"""
        variables = []
        for i in range(self.n_bits):
            variables.append(('x', i, f'x{i}'))
        for i in range(self.n_bits):
            variables.append(('y', i, f'y{i}'))

        monomials = []
        monomials_info = []

        # 常数项
        monomials.append([])
        monomials_info.append({
            'degree': 0, 'name': '1', 'type': 'constant',
            'x_bits': [], 'y_bits': [], 'index': 0
        })

        # 一次项
        for var in variables:
            monomials.append([var])
            info = {
                'degree': 1, 'name': var[2], 'type': 'linear',
                'x_bits': [var[1]] if var[0] == 'x' else [],
                'y_bits': [var[1]] if var[0] == 'y' else [],
                'index': len(monomials) - 1
            }
            monomials_info.append(info)

        # 二次项
        if max_degree >= 2:
            for i in range(len(variables)):
                for j in range(i + 1, len(variables)):
                    monomials.append([variables[i], variables[j]])
                    x_bits = []
                    y_bits = []
                    for var in [variables[i], variables[j]]:
                        if var[0] == 'x':
                            x_bits.append(var[1])
                        else:
                            y_bits.append(var[1])
                    info = {
                        'degree': 2,
                        'name': f"{variables[i][2]}{variables[j][2]}",
                        'type': 'quadratic',
                        'x_bits': x_bits,
                        'y_bits': y_bits,
                        'index': len(monomials) - 1
                    }
                    monomials_info.append(info)

        return monomials, monomials_info

    @staticmethod
    def evaluate_monomial(monomial, x_val, y_val):
        """计算单项式在(x,y)处的值（GF(2)上）"""
        if not monomial:
            return 1
        result = 1
        for var_type, bit_pos, _ in monomial:
            if var_type == 'x':
                result &= (x_val >> bit_pos) & 1
            else:
                result &= (y_val >> bit_pos) & 1
        return result

    def build_equation_system(self, truth_table, monomials):
        """构建线性方程组"""
        equations = []
        for x_val in range(self.n_inputs):
            y_val = truth_table[x_val]
            equations.append([self.evaluate_monomial(m, x_val, y_val) for m in monomials])
        return equations

    @staticmethod
    def gaussian_elimination_gf2(matrix):
        """GF(2)上的高斯消元法"""
        a_mat = [row[:] for row in matrix]
        n_rows = len(a_mat)
        n_cols = len(a_mat[0])

        pivot_rows = []
        pivot_cols = []

        current_row = 0
        for col in range(n_cols):
            pivot_found = False
            for i in range(current_row, n_rows):
                if a_mat[i][col] == 1:
                    a_mat[current_row], a_mat[i] = a_mat[i], a_mat[current_row]
                    pivot_found = True
                    pivot_rows.append(current_row)
                    pivot_cols.append(col)
                    break
            if pivot_found:
                for i in range(n_rows):
                    if i != current_row and a_mat[i][col] == 1:
                        for j in range(col, n_cols):
                            a_mat[i][j] ^= a_mat[current_row][j]
                current_row += 1
            if current_row >= n_rows:
                break

        return a_mat, pivot_rows, pivot_cols

    def find_null_space_gf2(self, matrix):
        """求解零空间"""
        rref, pivot_rows, pivot_cols = self.gaussian_elimination_gf2(matrix)
        n_cols = len(matrix[0])
        rank = len(pivot_cols)
        free_cols = [j for j in range(n_cols) if j not in pivot_cols]

        null_space = []
        for free_col in free_cols:
            solution = [0] * n_cols
            solution[free_col] = 1
            for i in range(rank - 1, -1, -1):
                pivot_row = pivot_rows[i]
                pivot_col = pivot_cols[i]
                value = 0
                for j in range(pivot_col + 1, n_cols):
                    if rref[pivot_row][j] == 1:
                        value ^= solution[j]
                solution[pivot_col] = value
            null_space.append(solution)

        return null_space, rank, free_cols

    def solve_implicit_equation(self, truth_table, max_degree=2):
        """求解隐式方程"""
        monomials, monomials_info = self.generate_all_monomials(max_degree)
        equations = self.build_equation_system(truth_table, monomials)
        null_space, rank, free_cols = self.find_null_space_gf2(equations)

        if null_space:
            non_zero = [s for s in null_space if any(c == 1 for c in s)]
            optimal = min(non_zero, key=lambda s: sum(s)) if non_zero else None
        else:
            optimal = None

        return optimal, null_space, monomials_info, monomials, rank

    @staticmethod
    def coefficients_to_expression(coefficients, monomials_info):
        """系数向量转方程字符串"""
        if not coefficients:
            return "0 = 0"
        terms = [info['name'] for c, info in zip(coefficients, monomials_info) if c == 1]
        return " ⊕ ".join(terms) + " = 0" if terms else "0 = 0"


class LinearAttackEngine:
    """线性攻击分析引擎"""

    def __init__(self, n_bits=3):
        self.n_bits = n_bits
        self.n_inputs = 1 << n_bits
        self.calculator = SBoxImplicitCalculator(n_bits)

    @staticmethod
    def extract_linear_approximations(solution_vector, monomials_info):
        """从解向量提取线性逼近"""
        x_bits_list = []
        y_bits_list = []
        has_constant = False

        for coeff, info in zip(solution_vector, monomials_info):
            if coeff == 1:
                if info['degree'] == 0:
                    has_constant = True
                elif info['degree'] == 1:
                    x_bits_list.extend(info['x_bits'])
                    y_bits_list.extend(info['y_bits'])

        x_bits_list = sorted(set(x_bits_list))
        y_bits_list = sorted(set(y_bits_list))

        if not x_bits_list and not y_bits_list:
            return []

        input_mask = sum(1 << b for b in x_bits_list)
        output_mask = sum(1 << b for b in y_bits_list)

        parts = []
        if has_constant:
            parts.append("1")
        parts.extend([f"x{b}" for b in x_bits_list])
        parts.extend([f"y{b}" for b in y_bits_list])
        expr = " ⊕ ".join(parts) + " = 0"

        return [{
            'input_mask': input_mask,
            'output_mask': output_mask,
            'x_bits': x_bits_list,
            'y_bits': y_bits_list,
            'has_constant': has_constant,
            'expression': expr
        }]

    def compute_bias(self, truth_table, input_mask, output_mask, has_constant=False):
        """计算线性偏差"""
        count = 0
        for x_val in range(self.n_inputs):
            y_val = truth_table[x_val]
            input_parity = bin(x_val & input_mask).count('1') % 2
            output_parity = bin(y_val & output_mask).count('1') % 2
            result = input_parity ^ output_parity
            expected = 1 if has_constant else 0
            if result == expected:
                count += 1

        half = self.n_inputs // 2
        bias = abs(count - half)
        probability = count / self.n_inputs
        return bias, probability, count

    def analyze_sbox(self, truth_table):
        """分析S盒的线性脆弱性"""
        optimal, null_space, monomials_info, _, _ = self.calculator.solve_implicit_equation(truth_table)

        if not optimal:
            return None

        all_approximations = []
        for solution in null_space:
            if not any(c == 1 for c in solution):
                continue
            approximations = self.extract_linear_approximations(solution, monomials_info)
            for approx in approximations:
                bias, prob, count = self.compute_bias(
                    truth_table, approx['input_mask'], approx['output_mask'], approx['has_constant']
                )
                approx['bias'] = bias
                approx['probability'] = prob
                approx['count'] = count
                all_approximations.append(approx)

        seen = set()
        unique = []
        for a in all_approximations:
            key = (a['input_mask'], a['output_mask'], a['has_constant'])
            if key not in seen:
                seen.add(key)
                unique.append(a)
        unique.sort(key=lambda x: x['bias'], reverse=True)

        best = unique[0] if unique else None
        max_bias = best['bias'] if best else 0
        half = self.n_inputs // 2

        if max_bias == half:
            level = "可破解(完美线性)"
        elif max_bias >= half * 3 // 4:
            level = "极易破解"
        elif max_bias >= half // 2:
            level = "可破解"
        elif max_bias > 0:
            level = "较难破解"
        else:
            level = "安全(无线性弱点)"

        return {
            'optimal_coefficients': optimal,
            'expression': self.calculator.coefficients_to_expression(optimal, monomials_info),
            'approximations': unique,
            'best': best,
            'max_bias': max_bias,
            'attack_level': level,
            'null_space_dim': len(null_space),
            'monomials_info': monomials_info
        }


# =====================================================================
# S盒性能指标计算函数
# =====================================================================

def calc_hamming_weight(val):
    """计算汉明重量"""
    count = 0
    while val:
        count += val & 1
        val >>= 1
    return count


def calc_inner_product(a, b):
    """计算GF(2)内积"""
    return calc_hamming_weight(a & b) % 2


def sbox_to_component_truth_tables(sbox, n_bits):
    """拆分S盒为各个分量的真值表"""
    size = 1 << n_bits
    components = []
    for bit in range(n_bits):
        tt = [(sbox[x] >> (n_bits - 1 - bit)) & 1 for x in range(size)]
        components.append(tt)
    return components


def mobius_transform(truth_table, n_bits):
    """Möbius变换求ANF系数"""
    size = 1 << n_bits
    anf = truth_table[:]
    for i in range(n_bits):
        step = 1 << i
        for j in range(size):
            if j & step:
                anf[j] ^= anf[j ^ step]
    return anf


def walsh_hadamard_transform(truth_table, n_bits):
    """Walsh-Hadamard变换"""
    size = 1 << n_bits
    wht = [1 - 2 * x for x in truth_table]
    for i in range(n_bits):
        step = 1 << i
        for j in range(size):
            if not (j & step):
                u_val, v_val = wht[j], wht[j | step]
                wht[j], wht[j | step] = u_val + v_val, u_val - v_val
    return wht


def nonlinearity_boolean(truth_table, n_bits):
    """布尔函数非线性度"""
    wht = walsh_hadamard_transform(truth_table, n_bits)
    return (1 << (n_bits - 1)) - max(abs(w) for w in wht) // 2


def nonlinearity_sbox(sbox, n_bits):
    """S盒非线性度"""
    size = 1 << n_bits
    min_nl = float('inf')
    for v in range(1, 1 << n_bits):
        tt = [calc_inner_product(v, sbox[x]) for x in range(size)]
        nl = nonlinearity_boolean(tt, n_bits)
        if nl < min_nl:
            min_nl = nl
    return int(min_nl)


def differential_uniformity(sbox, n_bits):
    """差分均匀度"""
    size = 1 << n_bits
    ddt = [[0] * size for _ in range(size)]
    for dx in range(size):
        for x in range(size):
            ddt[dx][sbox[x] ^ sbox[x ^ dx]] += 1
    return max(ddt[dx][dy] for dx in range(1, size) for dy in range(size))


def max_linear_approximation(sbox, n_bits):
    """最大线性逼近值"""
    size_in = 1 << n_bits
    lat = [[0] * size_in for _ in range(size_in)]
    for a in range(size_in):
        for b in range(size_in):
            count = sum(1 for x in range(size_in)
                        if calc_inner_product(a, x) == calc_inner_product(b, sbox[x]))
            lat[a][b] = count - (size_in >> 1)
    return max(abs(lat[a][b]) for a in range(1 << n_bits) for b in range(1, 1 << n_bits))


def algebraic_degree_sbox(sbox, n_bits):
    """S盒代数次数"""
    components = sbox_to_component_truth_tables(sbox, n_bits)
    max_deg = 0
    for comp in components:
        anf = mobius_transform(comp, n_bits)
        deg = max((calc_hamming_weight(k) for k in range(1 << n_bits) if anf[k] == 1), default=0)
        max_deg = max(max_deg, deg)
    return max_deg


def fixed_points_count(sbox, n_bits):
    """不动点个数"""
    return sum(1 for x in range(1 << n_bits) if sbox[x] == x)


def sac_property(sbox, n_bits):
    """严格雪崩准则"""
    size = 1 << n_bits
    components = sbox_to_component_truth_tables(sbox, n_bits)
    distances = []
    for tt in components:
        max_dist = 0
        for bit in range(n_bits):
            dist = abs(sum(1 for x in range(size) if tt[x] ^ tt[x ^ (1 << bit)]) - (size >> 1))
            max_dist = max(max_dist, dist)
        distances.append(max_dist)
    return distances


def is_involution(sbox, n_bits):
    """是否对合"""
    return all(sbox[sbox[x]] == x for x in range(1 << n_bits))


def compute_bit_independence(sbox, n_bits):
    """计算位独立性指标"""
    size = 1 << n_bits
    components = sbox_to_component_truth_tables(sbox, n_bits)
    correlations = []
    for i in range(n_bits):
        for j in range(i + 1, n_bits):
            tt1, tt2 = components[i], components[j]
            corr = sum(1 for x in range(size) if tt1[x] == tt2[x])
            correlations.append(abs(corr - size // 2))
    return max(correlations) if correlations else 0


def compute_avalanche_score(sbox, n_bits):
    """计算雪崩效应得分"""
    size = 1 << n_bits
    total_changes = 0
    total_pairs = 0
    for dx in range(1, size):
        for x in range(size):
            diff = sbox[x] ^ sbox[x ^ dx]
            total_changes += calc_hamming_weight(diff)
            total_pairs += 1
    return total_changes / total_pairs / n_bits


# =====================================================================
# 格式化函数
# =====================================================================

def format_percentage(count, total):
    """格式化百分比"""
    if total == 0:
        return "100.00%"
    return f"{count / total * 100:.2f}%"


# =====================================================================
# 3×3 S盒处理与优选
# =====================================================================

def process_3bit_sboxes():
    """处理所有40320个3-bit双射S盒并进行优选"""
    n_bits = 3
    print("\n" + "=" * 80)
    print(f"处理{n_bits}x{n_bits} S盒（共40320个）- 含优选分析")
    print("=" * 80)

    calculator = SBoxImplicitCalculator(n_bits=n_bits)
    engine = LinearAttackEngine(n_bits=n_bits)

    all_results = []
    start_time = time.time()

    all_perms = list(permutations(range(8)))
    total = len(all_perms)

    for sbox_id, perm in enumerate(all_perms, 1):
        truth_table = list(perm)
        optimal, null_space, monomials_info, _, _ = calculator.solve_implicit_equation(truth_table)
        expression = calculator.coefficients_to_expression(optimal, monomials_info) if optimal else "无解"

        nl = nonlinearity_sbox(truth_table, n_bits)
        diff = differential_uniformity(truth_table, n_bits)
        deg = algebraic_degree_sbox(truth_table, n_bits)
        fp = fixed_points_count(truth_table, n_bits)
        max_lat = max_linear_approximation(truth_table, n_bits)
        sac = sac_property(truth_table, n_bits)
        inv = is_involution(truth_table, n_bits)
        bit_ind = compute_bit_independence(truth_table, n_bits)
        ava_score = compute_avalanche_score(truth_table, n_bits)

        analysis = engine.analyze_sbox(truth_table)

        all_results.append({
            'id': sbox_id,
            'truth_table': truth_table,
            'expression': expression,
            'optimal_coefficients': optimal,
            'null_space': null_space,
            'monomials_info': monomials_info,
            'nl': nl,
            'diff': diff,
            'deg': deg,
            'fp': fp,
            'max_lat': max_lat,
            'sac': sac,
            'inv': inv,
            'bit_independence': bit_ind,
            'avalanche_score': ava_score,
            'attack_level': analysis['attack_level'] if analysis else 'N/A'
        })

        if sbox_id % 5000 == 0:
            elapsed = time.time() - start_time
            print(f"  进度: {sbox_id:6d}/{total} ({format_percentage(sbox_id, total)}) - {elapsed:.1f}秒")

    print(f"  完成！总耗时: {time.time() - start_time:.1f}秒")

    # 解向量分类
    solution_classes = defaultdict(list)
    for r in all_results:
        if r['optimal_coefficients']:
            solution_classes[tuple(r['optimal_coefficients'])].append(r)

    print(f"\n  解向量类别总数: {len(solution_classes)}")

    # ===== 优选分析 =====
    print("\n  进行优选分析...")

    # 定义评分标准
    def calculate_score(r):
        """计算S盒综合评分"""
        score = 0

        # NL: 2分满分，越低越好（2为最优）
        nl_score = {2: 10, 1: 4, 0: 1}
        score += nl_score.get(r['nl'], 0)

        # δ: 2分满分，越低越好（2为最优）
        diff_score = {2: 10, 4: 4, 6: 1, 8: 0}
        score += diff_score.get(r['diff'], 0)

        # deg: 2分满分，越高越好
        deg_score = {2: 8, 1: 3, 0: 0}
        score += deg_score.get(r['deg'], 0)

        # |LAT_max|: 越低越好（0为最优）
        lat_score = {0: 8, 2: 5, 4: 2}
        score += lat_score.get(r['max_lat'], 0)

        # 不动点: 越少越好
        fp_score = {0: 5, 2: 3, 4: 1, 8: 0}
        score += fp_score.get(r['fp'], 0)

        # 安全等级
        attack_score = {
            '安全(无线性弱点)': 10,
            '较难破解': 6,
            '可破解': 3,
            '极易破解': 1,
            '可破解(完美线性)': 0
        }
        score += attack_score.get(r['attack_level'], 0)

        # 位独立性: 越接近0越好
        if r['bit_independence'] <= 1:
            score += 4
        elif r['bit_independence'] <= 2:
            score += 2

        # 雪崩效应: 越接近0.5越好
        if 0.45 <= r['avalanche_score'] <= 0.55:
            score += 4
        elif 0.40 <= r['avalanche_score'] <= 0.60:
            score += 2

        return score

    # 对每个解向量类别进行评分
    class_scores = []
    for coeff_key, sboxes in solution_classes.items():
        # 计算类内平均分
        scores = [calculate_score(r) for r in sboxes]
        avg_score = sum(scores) / len(scores)
        max_score = max(scores)
        min_score = min(scores)

        # 类内指标统计
        nls = [r['nl'] for r in sboxes]
        diffs = [r['diff'] for r in sboxes]
        degs = [r['deg'] for r in sboxes]
        lats = [r['max_lat'] for r in sboxes]
        fps = [r['fp'] for r in sboxes]

        # 类内最优S盒
        best_sbox = max(sboxes, key=lambda r: calculate_score(r))

        # 安全性统计
        secure_count = sum(1 for r in sboxes if r['attack_level'] == '安全(无线性弱点)')

        class_scores.append({
            'coeff_key': coeff_key,
            'num_sboxes': len(sboxes),
            'avg_score': avg_score,
            'max_score': max_score,
            'min_score': min_score,
            'best_sbox': best_sbox,
            'nl_range': (min(nls), max(nls)),
            'diff_range': (min(diffs), max(diffs)),
            'deg_range': (min(degs), max(degs)),
            'lat_range': (min(lats), max(lats)),
            'fp_range': (min(fps), max(fps)),
            'secure_count': secure_count,
            'secure_ratio': secure_count / len(sboxes),
            'example': sboxes[0]
        })

    # 按平均分排序
    class_scores.sort(key=lambda x: x['avg_score'], reverse=True)

    # 保存隐式方程
    print("\n  保存 3x3.txt...")
    with open('3x3.txt', 'w', encoding='utf-8') as f:
        f.write("3-bit双射S盒隐式方程\n")
        f.write(f"总数: {total}\n")
        f.write("形式: F(x0,x1,x2,y0,y1,y2) = 0, 次数<3\n")
        f.write("=" * 80 + "\n\n")
        for r in all_results:
            f.write(f"S盒 #{r['id']:5d}: {r['expression']}\n")

    # 保存优选结果
    print("  保存优选结果到 3x3_优选.txt...")

    with open('3x3_优选.txt', 'w', encoding='utf-8') as f:
        f.write("3-bit双射S盒优选分析报告\n")
        f.write("从解向量分类中基于密码学性能指标进一步优选\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"总S盒数量: {total}\n")
        f.write(f"解向量类别数: {len(solution_classes)}\n")
        f.write(f"分析日期: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # 评分标准说明
        f.write("评分标准（满分59分）:\n")
        f.write("-" * 40 + "\n")
        f.write("  非线性度(NL): NL=2得10分, NL=1得4分, NL=0得1分\n")
        f.write("  差分均匀度(δ): δ=2得10分, δ=4得4分, δ=6得1分, δ=8得0分\n")
        f.write("  代数次数(deg): deg=2得8分, deg=1得3分, deg=0得0分\n")
        f.write("  |LAT_max|: =0得8分, =2得5分, =4得2分\n")
        f.write("  不动点: =0得5分, =2得3分, =4得1分, =8得0分\n")
        f.write("  安全等级: 安全得10分, 较难得6分, 可破解得3分, 极易得1分, 完美线性得0分\n")
        f.write("  位独立性: ≤1得4分, ≤2得2分\n")
        f.write("  雪崩效应: [0.45,0.55]得4分, [0.40,0.60]得2分\n\n")

        # 分数段分布
        f.write("各类评分分布:\n")
        f.write("-" * 40 + "\n")
        score_ranges = defaultdict(int)
        for cs in class_scores:
            s = cs['avg_score']
            if s >= 50:
                score_ranges['50-59 (优秀)'] += 1
            elif s >= 40:
                score_ranges['40-49 (良好)'] += 1
            elif s >= 30:
                score_ranges['30-39 (中等)'] += 1
            elif s >= 20:
                score_ranges['20-29 (较差)'] += 1
            else:
                score_ranges['<20 (差)'] += 1

        for range_name in ['50-59 (优秀)', '40-49 (良好)', '30-39 (中等)', '20-29 (较差)', '<20 (差)']:
            count = score_ranges.get(range_name, 0)
            if count > 0:
                f.write(f"  {range_name}: {count}类 ({format_percentage(count, len(class_scores))})\n")

        f.write("\n" + "=" * 80 + "\n\n")

        # 前50类详细分析
        f.write("优选排名（前50类）:\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"{'排名':>4} {'数量':>6} {'占比':>8} {'评分':>6} {'安全%':>8} "
                f"{'NL':>8} {'δ':>8} {'deg':>6} {'|LAT|':>6} {'FP':>6} {'雪崩':>6}\n")
        f.write("-" * 85 + "\n")

        for rank, cs in enumerate(class_scores[:50], 1):
            pct = format_percentage(cs['num_sboxes'], total)
            secure_pct = format_percentage(cs['secure_count'], cs['num_sboxes'])
            nl_str = f"{cs['nl_range'][0]}-{cs['nl_range'][1]}"
            diff_str = f"{cs['diff_range'][0]}-{cs['diff_range'][1]}"
            deg_str = f"{cs['deg_range'][0]}-{cs['deg_range'][1]}"
            lat_str = f"{cs['lat_range'][0]}-{cs['lat_range'][1]}"
            fp_str = f"{cs['fp_range'][0]}-{cs['fp_range'][1]}"

            f.write(f"{rank:>4} {cs['num_sboxes']:>6} {pct:>8} {cs['avg_score']:>6.1f} "
                    f"{secure_pct:>8} {nl_str:>8} {diff_str:>8} {deg_str:>6} "
                    f"{lat_str:>6} {fp_str:>6} {cs['best_sbox']['avalanche_score']:>6.3f}\n")

        # 前10类详细展示
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("前10类详细展示\n")
        f.write("=" * 80 + "\n\n")

        for rank, cs in enumerate(class_scores[:10], 1):
            f.write(f"排名 {rank} (评分: {cs['avg_score']:.1f}/59)\n")
            f.write("=" * 50 + "\n")

            example = cs['example']
            best = cs['best_sbox']

            f.write(f"类别大小: {cs['num_sboxes']}个 ({format_percentage(cs['num_sboxes'], total)})\n")
            f.write(f"解向量: [{', '.join(str(c) for c in cs['coeff_key'])}]\n")
            f.write(f"隐式方程: {example['expression']}\n\n")

            f.write(f"性能指标:\n")
            f.write(f"  NL范围: {cs['nl_range']}\n")
            f.write(f"  δ范围: {cs['diff_range']}\n")
            f.write(f"  deg范围: {cs['deg_range']}\n")
            f.write(f"  |LAT|范围: {cs['lat_range']}\n")
            f.write(f"  FP范围: {cs['fp_range']}\n")
            f.write(
                f"  安全比例: {cs['secure_count']}/{cs['num_sboxes']} ({format_percentage(cs['secure_count'], cs['num_sboxes'])})\n\n")

            f.write(f"类内最优S盒 (#{best['id']}):\n")
            f.write(f"  真值表: {best['truth_table']}\n")
            f.write(f"  NL={best['nl']}, δ={best['diff']}, deg={best['deg']}\n")
            f.write(f"  |LAT_max|={best['max_lat']}, FP={best['fp']}\n")
            f.write(f"  位独立性={best['bit_independence']}, 雪崩={best['avalanche_score']:.3f}\n")
            f.write(f"  攻击等级: {best['attack_level']}\n")
            f.write(f"  评分: {calculate_score(best)}/59\n\n")

            f.write("-" * 50 + "\n\n")

        # 最优S盒推荐
        f.write("\n" + "=" * 80 + "\n")
        f.write("最优S盒推荐\n")
        f.write("=" * 80 + "\n\n")

        # 找出所有满分或接近满分的S盒
        top_sboxes = []
        for r in all_results:
            score = calculate_score(r)
            if score >= 50:
                top_sboxes.append((score, r))

        top_sboxes.sort(key=lambda x: x[0], reverse=True)

        f.write(f"评分≥50的S盒: {len(top_sboxes)}个\n\n")

        for rank, (score, r) in enumerate(top_sboxes[:20], 1):
            f.write(f"推荐 #{rank} (评分: {score}/59)\n")
            f.write(f"  S盒 ID: {r['id']}\n")
            f.write(f"  真值表: {r['truth_table']}\n")
            f.write(f"  方程: {r['expression']}\n")
            f.write(f"  NL={r['nl']}, δ={r['diff']}, deg={r['deg']}\n")
            f.write(f"  |LAT|={r['max_lat']}, FP={r['fp']}\n")
            f.write(f"  位独立性={r['bit_independence']}, 雪崩={r['avalanche_score']:.3f}\n")
            f.write(f"  攻击等级: {r['attack_level']}\n\n")

        # 安全级别统计
        f.write("\n安全级别分布:\n")
        f.write("-" * 40 + "\n")
        attack_dist = defaultdict(int)
        for r in all_results:
            attack_dist[r['attack_level']] += 1
        for level in ['安全(无线性弱点)', '较难破解', '可破解', '极易破解', '可破解(完美线性)']:
            count = attack_dist.get(level, 0)
            if count > 0:
                f.write(f"  {level}: {count}个 ({format_percentage(count, total)})\n")

    print("  3x3 S盒处理完成！")
    return all_results, class_scores


# =====================================================================
# 4×4 S盒处理
# =====================================================================

def get_4bit_affine_representatives():
    """4×4 S盒的16个仿射等价类代表元"""
    return [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        [12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2],
        [0, 1, 2, 3, 4, 6, 7, 5, 8, 9, 10, 11, 12, 13, 15, 14],
        [0, 1, 2, 4, 3, 5, 6, 7, 8, 9, 10, 12, 11, 13, 14, 15],
        [0, 2, 1, 3, 4, 5, 6, 7, 8, 10, 9, 11, 12, 13, 14, 15],
        [0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15],
        [0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15],
        [0, 2, 4, 6, 1, 3, 5, 7, 8, 10, 12, 14, 9, 11, 13, 15],
        [0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11, 12, 13, 14, 15],
        [0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11, 12, 13, 14, 15],
        [0, 1, 2, 4, 3, 5, 8, 9, 6, 7, 10, 12, 11, 13, 14, 15],
        [0, 1, 4, 5, 2, 3, 8, 9, 6, 7, 12, 13, 10, 11, 14, 15],
        [0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 15],
        [0, 1, 2, 4, 3, 6, 5, 8, 7, 9, 10, 12, 11, 14, 13, 15],
        [0, 2, 4, 6, 1, 3, 5, 8, 7, 9, 11, 13, 10, 12, 14, 15],
        [0, 2, 1, 4, 3, 6, 8, 10, 7, 9, 5, 12, 11, 14, 13, 15],
    ]


def process_4bit_sboxes():
    """处理16个4-bit仿射等价类代表元"""
    n_bits = 4
    total = 16

    print("\n" + "=" * 80)
    print(f"处理{n_bits}x{n_bits} S盒（16个仿射等价类代表元）")
    print("=" * 80)

    calculator = SBoxImplicitCalculator(n_bits=n_bits)
    engine = LinearAttackEngine(n_bits=n_bits)

    representatives = get_4bit_affine_representatives()
    all_results = []

    known_sboxes = {
        'PRESENT': [0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2],
        'GIFT': [0x1, 0xA, 0x4, 0xC, 0x6, 0xF, 0x3, 0x9, 0x2, 0xD, 0xB, 0x7, 0x5, 0x0, 0x8, 0xE],
        'Midori': [0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6]
    }

    for sbox_id, truth_table in enumerate(representatives, 1):
        optimal, null_space, monomials_info, monomials, rank = calculator.solve_implicit_equation(truth_table)
        expression = calculator.coefficients_to_expression(optimal, monomials_info) if optimal else "无解"

        nl = nonlinearity_sbox(truth_table, n_bits)
        diff = differential_uniformity(truth_table, n_bits)
        deg = algebraic_degree_sbox(truth_table, n_bits)
        fp = fixed_points_count(truth_table, n_bits)
        max_lat = max_linear_approximation(truth_table, n_bits)
        sac = sac_property(truth_table, n_bits)
        inv = is_involution(truth_table, n_bits)
        bit_ind = compute_bit_independence(truth_table, n_bits)
        ava_score = compute_avalanche_score(truth_table, n_bits)

        analysis = engine.analyze_sbox(truth_table)

        all_results.append({
            'id': sbox_id,
            'truth_table': truth_table,
            'expression': expression,
            'optimal_coefficients': optimal,
            'num_monomials': len(monomials),
            'rank': rank,
            'nullity': len(null_space),
            'nl': nl,
            'diff': diff,
            'deg': deg,
            'fp': fp,
            'max_lat': max_lat,
            'sac': sac,
            'inv': inv,
            'bit_independence': bit_ind,
            'avalanche_score': ava_score,
            'attack_level': analysis['attack_level'] if analysis else 'N/A'
        })

        print(f"  S盒 #{sbox_id}: NL={nl}, δ={diff}, deg={deg}")

    # 保存隐式方程
    print("\n  保存 4x4.txt...")
    with open('4x4.txt', 'w', encoding='utf-8') as f:
        f.write("4-bit仿射等价类代表元隐式方程\n")
        f.write(f"总数: {total}个仿射等价类\n")
        f.write("形式: F(x0,x1,x2,x3,y0,y1,y2,y3) = 0, 次数<3\n")
        f.write("=" * 80 + "\n\n")
        for r in all_results:
            f.write(f"等价类 #{r['id']:2d}:\n")
            f.write(f"  真值表: {r['truth_table']}\n")
            f.write(f"  方程: {r['expression']}\n\n")

    # 保存统计报告
    print("  保存统计报告到 4x4_统计.txt...")

    with open('4x4_统计.txt', 'w', encoding='utf-8') as f:
        f.write("4-bit仿射等价类代表元性能统计报告\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"总数: {total}个仿射等价类\n\n")

        f.write("1. 各等价类完整性能指标\n")
        f.write("-" * 75 + "\n")
        f.write(f"{'ID':>3} {'NL':>4} {'δ':>4} {'deg':>4} {'FP':>4} {'|LAT|':>6} "
                f"{'位独立':>6} {'雪崩':>6} {'对合':>4} {'攻击等级':<20}\n")
        f.write("-" * 75 + "\n")
        for r in all_results:
            f.write(f"{r['id']:>3} {r['nl']:>4} {r['diff']:>4} {r['deg']:>4} "
                    f"{r['fp']:>4} {r['max_lat']:>6} {r['bit_independence']:>6} "
                    f"{r['avalanche_score']:>6.3f} {'是' if r['inv'] else '否':>4} "
                    f"{r['attack_level']:<20}\n")

        f.write("\n2. 已知密码S盒对比\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'S盒':>12} {'NL':>4} {'δ':>4} {'deg':>4} {'FP':>4} {'|LAT|':>6} {'对合':>4}\n")
        f.write("-" * 50 + "\n")
        for name, sbox in known_sboxes.items():
            nl_val = nonlinearity_sbox(sbox, n_bits)
            diff_val = differential_uniformity(sbox, n_bits)
            deg_val = algebraic_degree_sbox(sbox, n_bits)
            fp_val = fixed_points_count(sbox, n_bits)
            lat_val = max_linear_approximation(sbox, n_bits)
            inv_val = is_involution(sbox, n_bits)
            f.write(f"{name:>12} {nl_val:>4} {diff_val:>4} {deg_val:>4} "
                    f"{fp_val:>4} {lat_val:>6} {'是' if inv_val else '否':>4}\n")

        f.write("\n3. 统计汇总\n")
        f.write("-" * 40 + "\n")
        nl_vals = [r['nl'] for r in all_results]
        diff_vals = [r['diff'] for r in all_results]
        deg_vals = [r['deg'] for r in all_results]
        lat_vals = [r['max_lat'] for r in all_results]
        fp_vals = [r['fp'] for r in all_results]

        f.write(f"{'指标':>15} {'min':>6} {'max':>6} {'avg':>8}\n")
        f.write("-" * 35 + "\n")
        f.write(f"{'非线性度':>15} {min(nl_vals):>6} {max(nl_vals):>6} {sum(nl_vals) / total:>8.2f}\n")
        f.write(f"{'差分均匀度':>15} {min(diff_vals):>6} {max(diff_vals):>6} {sum(diff_vals) / total:>8.2f}\n")
        f.write(f"{'代数次数':>15} {min(deg_vals):>6} {max(deg_vals):>6} {sum(deg_vals) / total:>8.2f}\n")
        f.write(f"{'|LAT_max|':>15} {min(lat_vals):>6} {max(lat_vals):>6} {sum(lat_vals) / total:>8.2f}\n")
        f.write(f"{'不动点':>15} {min(fp_vals):>6} {max(fp_vals):>6} {sum(fp_vals) / total:>8.2f}\n")

        inv_count = sum(1 for r in all_results if r['inv'])
        no_fp_count = sum(1 for r in all_results if r['fp'] == 0)
        f.write(f"\n对合等价类: {inv_count}个 ({format_percentage(inv_count, total)})\n")
        f.write(f"无不动点等价类: {no_fp_count}个 ({format_percentage(no_fp_count, total)})\n")

    print("  4x4 S盒处理完成！")
    return all_results


# =====================================================================
# 主程序
# =====================================================================

def main():
    """主程序"""
    print("=" * 80)
    print("S盒隐式方程计算与密码学性能优选系统")
    print("支持3x3和4x4 S盒")
    print("=" * 80)

    start_time = time.time()

    results_3, class_scores_3 = process_3bit_sboxes()
    results_4 = process_4bit_sboxes()

    total_time = time.time() - start_time

    print("\n" + "=" * 80)
    print(f"全部处理完成！总耗时: {total_time:.1f}秒")
    print("=" * 80)
    print("\n生成的文件:")
    print("  3x3 S盒:")
    print("    3x3.txt       - 所有40320个S盒的隐式方程")
    print("    3x3_优选.txt  - 解向量分类优选分析（含评分排名）")
    print("  4x4 S盒:")
    print("    4x4.txt       - 16个仿射等价类的隐式方程")
    print("    4x4_统计.txt  - 性能统计报告")

    print(f"\n统计摘要:")
    print(f"  3x3: {len(results_3)}个S盒, {len(class_scores_3)}个解向量类别")
    print(f"  4x4: {len(results_4)}个仿射等价类")

    # 打印优选摘要
    if class_scores_3:
        print(f"\n3x3优选摘要:")
        print(f"  最高评分类: {class_scores_3[0]['avg_score']:.1f}分")
        print(f"  前10类平均分: {sum(c['avg_score'] for c in class_scores_3[:10]) / 10:.1f}分")



# =====================================================================
# 论文补充实验模块：ANF / DDT / LAT / 线性攻击演示
# 说明：
#   1. 不替代原有统计程序，只补充可直接写入论文的实验输出。
#   2. 建议运行：python main_enhanced.py --supplement
#   3. 如需保留原完整枚举，运行：python main_enhanced.py --full
# =====================================================================

import argparse


def parity(value):
    """返回整数二进制表示中1的个数模2"""
    return bin(value).count("1") & 1


def bit_name(bit_index, n_bits, prefix="x"):
    """统一变量命名：bit_index=0 表示最低位；论文中可根据需要改成x1,x2形式"""
    return f"{prefix}{bit_index}"


def anf_to_expression(anf_coeffs, n_bits, var_prefix="x"):
    """将ANF系数转换成字符串表达式"""
    terms = []
    for mask, coeff in enumerate(anf_coeffs):
        if coeff == 0:
            continue
        if mask == 0:
            terms.append("1")
            continue
        factors = [bit_name(i, n_bits, var_prefix) for i in range(n_bits) if (mask >> i) & 1]
        terms.append("".join(factors))
    return " ⊕ ".join(terms) if terms else "0"


def component_anf_report(sbox, n_bits):
    """返回S盒每个分量函数的ANF表达式和次数"""
    components = sbox_to_component_truth_tables(sbox, n_bits)
    rows = []
    for idx, tt in enumerate(components):
        anf = mobius_transform(tt, n_bits)
        expr = anf_to_expression(anf, n_bits, "x")
        degree = max((calc_hamming_weight(mask) for mask, c in enumerate(anf) if c), default=0)
        rows.append({
            "component": f"y{idx}",
            "truth_table": tt,
            "anf_coeffs": anf,
            "anf_expr": expr,
            "degree": degree
        })
    return rows


def build_ddt(sbox, n_bits):
    """构造差分分布表DDT"""
    size = 1 << n_bits
    ddt = [[0] * size for _ in range(size)]
    for dx in range(size):
        for x in range(size):
            dy = sbox[x] ^ sbox[x ^ dx]
            ddt[dx][dy] += 1
    return ddt


def build_lat(sbox, n_bits):
    """构造线性分布表LAT，表项为 count - 2^(n-1)"""
    size = 1 << n_bits
    lat = [[0] * size for _ in range(size)]
    for a in range(size):
        for b in range(size):
            count = 0
            for x in range(size):
                if calc_inner_product(a, x) == calc_inner_product(b, sbox[x]):
                    count += 1
            lat[a][b] = count - (size >> 1)
    return lat


def mask_to_linear_expr(mask, n_bits, prefix):
    """把掩码转换成线性表达式"""
    terms = [f"{prefix}{i}" for i in range(n_bits) if (mask >> i) & 1]
    return " ⊕ ".join(terms) if terms else "0"


def best_linear_approximations_from_lat(lat, n_bits, top_k=10):
    """从LAT中选出最大线性逼近。排除a=b=0，并优先保留非零输入、非零输出掩码。"""
    size = 1 << n_bits
    items = []
    for a in range(size):
        for b in range(size):
            if a == 0 and b == 0:
                continue
            if a == 0 or b == 0:
                # 这类只反映平衡性，不作为线性攻击主逼近
                continue
            val = lat[a][b]
            if val == 0:
                continue
            prob = (val + (size >> 1)) / size
            bias = abs(val) / size
            items.append({
                "input_mask": a,
                "output_mask": b,
                "lat_value": val,
                "abs_lat": abs(val),
                "probability": prob,
                "bias": bias,
                "expr": f"{mask_to_linear_expr(a, n_bits, 'x')} = {mask_to_linear_expr(b, n_bits, 'y')}"
            })
    items.sort(key=lambda r: (r["abs_lat"], abs(r["probability"] - 0.5)), reverse=True)
    return items[:top_k]


def write_matrix(f, matrix, n_bits, title):
    """写矩阵表，适合DDT/LAT"""
    size = 1 << n_bits
    f.write(title + "\n")
    f.write("-" * 80 + "\n")
    header = "      " + " ".join(f"{i:X}".rjust(4) for i in range(size))
    f.write(header + "\n")
    for i, row in enumerate(matrix):
        f.write(f"{i:X}".rjust(4) + "  " + " ".join(str(v).rjust(4) for v in row) + "\n")
    f.write("\n")


def write_sbox_detailed_report(name, sbox, n_bits, out_dir):
    """为单个S盒输出可入论文的ANF、DDT、LAT、最佳线性逼近"""
    out_path = Path(out_dir) / f"{name}_{n_bits}bit_detailed_report.txt"
    ddt = build_ddt(sbox, n_bits)
    lat = build_lat(sbox, n_bits)
    anfs = component_anf_report(sbox, n_bits)
    best_las = best_linear_approximations_from_lat(lat, n_bits, top_k=12)

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"{name} S盒详细实验报告\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"真值表: {sbox}\n")
        f.write(f"规模: {n_bits}x{n_bits}\n\n")

        f.write("1. 核心密码学指标\n")
        f.write("-" * 80 + "\n")
        f.write(f"非线性度 NL = {nonlinearity_sbox(sbox, n_bits)}\n")
        f.write(f"差分均匀性 δ = {differential_uniformity(sbox, n_bits)}\n")
        f.write(f"代数次数 deg = {algebraic_degree_sbox(sbox, n_bits)}\n")
        f.write(f"最大线性表绝对值 |LAT_max| = {max_linear_approximation(sbox, n_bits)}\n")
        f.write(f"不动点个数 FP = {fixed_points_count(sbox, n_bits)}\n")
        f.write(f"是否对合 = {'是' if is_involution(sbox, n_bits) else '否'}\n")
        f.write(f"位独立性指标 = {compute_bit_independence(sbox, n_bits)}\n")
        f.write(f"平均雪崩效应 = {compute_avalanche_score(sbox, n_bits):.4f}\n\n")

        f.write("2. 各分量布尔函数ANF\n")
        f.write("-" * 80 + "\n")
        for row in anfs:
            f.write(f"{row['component']} = {row['anf_expr']}    deg={row['degree']}\n")
        f.write("\n")

        write_matrix(f, ddt, n_bits, "3. 差分分布表 DDT")
        write_matrix(f, lat, n_bits, "4. 线性分布表 LAT")

        f.write("5. LAT给出的主要线性逼近\n")
        f.write("-" * 80 + "\n")
        f.write("说明：LAT值为 count - 2^(n-1)，概率 = count / 2^n，偏差 = |P-1/2|。\n")
        for i, item in enumerate(best_las, 1):
            count = item["lat_value"] + (1 << (n_bits - 1))
            f.write(
                f"{i:02d}. a={item['input_mask']:X}, b={item['output_mask']:X}, "
                f"LAT={item['lat_value']:+d}, count={count}, "
                f"P={item['probability']:.4f}, bias={item['bias']:.4f}, "
                f"{item['expr']}\n"
            )
        f.write("\n")
    return out_path


def one_round_key_recovery_demo(sbox, n_bits, true_key, chosen_approx=None):
    """
    简化的一轮线性攻击演示：
    加密结构：C = S(P xor K)
    已知明文P、密文C，枚举K'，统计线性逼近
        <a, P xor K'> = <b, C>
    对所有K'的偏差排序，观察真实密钥是否排第一。

    该实验不是完整SPN攻击，但能说明LAT偏差如何转化为密钥区分器。
    """
    size = 1 << n_bits
    lat = build_lat(sbox, n_bits)
    best = chosen_approx or best_linear_approximations_from_lat(lat, n_bits, top_k=1)[0]
    a = best["input_mask"]
    b = best["output_mask"]

    pairs = []
    for p in range(size):
        c = sbox[p ^ true_key]
        pairs.append((p, c))

    results = []
    for guess_key in range(size):
        count = 0
        for p, c in pairs:
            left = calc_inner_product(a, p ^ guess_key)
            right = calc_inner_product(b, c)
            if left == right:
                count += 1
        bias = abs(count - (size >> 1))
        results.append({
            "guess_key": guess_key,
            "count": count,
            "bias": bias,
            "is_true_key": guess_key == true_key
        })

    results.sort(key=lambda r: r["bias"], reverse=True)
    rank = next(i + 1 for i, r in enumerate(results) if r["is_true_key"])
    return {
        "approx": best,
        "true_key": true_key,
        "rank": rank,
        "results": results
    }


def write_linear_attack_demo_report(name, sbox, n_bits, out_dir):
    """输出简化线性攻击实验报告"""
    out_path = Path(out_dir) / f"{name}_{n_bits}bit_linear_attack_demo.txt"
    true_keys = list(range(1 << n_bits))
    demos = [one_round_key_recovery_demo(sbox, n_bits, k) for k in true_keys]
    success_top1 = sum(1 for d in demos if d["rank"] == 1)
    success_top3 = sum(1 for d in demos if d["rank"] <= 3)

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"{name} S盒简化线性攻击实验\n")
        f.write("=" * 80 + "\n\n")
        f.write("实验结构: C = S(P xor K)\n")
        f.write("攻击方法: 根据LAT中最大偏差的线性逼近，枚举K'并统计偏差。\n")
        f.write("说明: 这是线性攻击原理演示，用于论文中“线性攻击应用”部分；完整多轮SPN攻击可作为后续扩展。\n\n")

        first = demos[0]
        approx = first["approx"]
        f.write("采用的线性逼近:\n")
        f.write(f"  a={approx['input_mask']:X}, b={approx['output_mask']:X}\n")
        f.write(f"  表达式: {approx['expr']}\n")
        f.write(f"  LAT={approx['lat_value']:+d}, P={approx['probability']:.4f}, bias={approx['bias']:.4f}\n\n")

        f.write("攻击成功率统计:\n")
        f.write("-" * 80 + "\n")
        f.write(f"密钥总数: {len(true_keys)}\n")
        f.write(f"真实密钥排名第1的次数: {success_top1}/{len(true_keys)} = {success_top1 / len(true_keys):.2%}\n")
        f.write(f"真实密钥排名前3的次数: {success_top3}/{len(true_keys)} = {success_top3 / len(true_keys):.2%}\n\n")

        f.write("逐密钥实验结果:\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'真实K':>6} {'排名':>6} {'第一候选':>8} {'第一偏差':>8} {'真实偏差':>8}\n")
        for demo in demos:
            top = demo["results"][0]
            true_row = next(r for r in demo["results"] if r["is_true_key"])
            f.write(
                f"{demo['true_key']:>6X} {demo['rank']:>6} "
                f"{top['guess_key']:>8X} {top['bias']:>8} {true_row['bias']:>8}\n"
            )
    return out_path


def run_supplemental_reports():
    """运行论文补充实验输出"""
    out_dir = Path("supplement_reports")
    out_dir.mkdir(exist_ok=True)

    sboxes_3bit = {
        "S3_optimal_example": [0, 1, 3, 6, 7, 4, 5, 2],
        "S3_weak_example": [5, 3, 6, 1, 0, 7, 2, 4],
    }

    sboxes_4bit = {
        "PRESENT": [0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2],
        "GIFT": [0x1, 0xA, 0x4, 0xC, 0x6, 0xF, 0x3, 0x9, 0x2, 0xD, 0xB, 0x7, 0x5, 0x0, 0x8, 0xE],
        "Midori": [0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6],
    }

    generated = []

    for name, sbox in sboxes_3bit.items():
        generated.append(write_sbox_detailed_report(name, sbox, 3, out_dir))
        generated.append(write_linear_attack_demo_report(name, sbox, 3, out_dir))

    for name, sbox in sboxes_4bit.items():
        generated.append(write_sbox_detailed_report(name, sbox, 4, out_dir))
        generated.append(write_linear_attack_demo_report(name, sbox, 4, out_dir))

    summary_path = out_dir / "supplement_experiment_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("补充实验建议与生成文件说明\n")
        f.write("=" * 80 + "\n\n")
        f.write("建议补充到论文中的实验:\n")
        f.write("1. ANF详细表：3bit选一个最优样例和一个弱样例，4bit选PRESENT、GIFT、Midori。\n")
        f.write("2. DDT完整表：3bit可全文展示，4bit建议正文展示PRESENT或Midori，其他放附录。\n")
        f.write("3. LAT完整表：用于支撑线性逼近分析，正文展示最大偏差行或典型S盒LAT。\n")
        f.write("4. 简化线性攻击实验：用 C=S(P xor K) 展示LAT偏差如何用于密钥区分。\n")
        f.write("5. 综合对比表：NL、δ、deg、|LAT_max|、FP、对合性、雪崩分数。\n\n")
        f.write("生成文件:\n")
        for p in generated:
            f.write(f"- {p.name}\n")
    generated.append(summary_path)

    print("\n补充实验报告已生成在 supplement_reports/ 目录下：")
    for p in generated:
        print(f"  {p}")


def cli_main():
    parser = argparse.ArgumentParser(description="S盒隐式方程与论文补充实验系统")
    parser.add_argument("--full", action="store_true", help="运行原完整实验：3bit全枚举 + 4bit代表元")
    parser.add_argument("--supplement", action="store_true", help="只运行论文补充实验：ANF/DDT/LAT/线性攻击演示")
    args = parser.parse_args()

    if args.full:
        main()
    elif args.supplement:
        run_supplemental_reports()
    else:
        print("未指定运行模式。建议：")
        print("  python main_enhanced.py --supplement   # 生成论文补充实验")
        print("  python main_enhanced.py --full         # 运行原完整枚举统计")


if __name__ == "__main__":
    cli_main()
