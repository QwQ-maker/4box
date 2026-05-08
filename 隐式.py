# ============================================================
# 新增模块A：隐式代数表示（Implicit Algebraic Representation）
# ============================================================

def compute_implicit_representation(sbox, n):
    """
    从ANF计算S盒的隐式代数表示。
    对每个分量方程 yi = fi(x1,...,xn)，移项得 fi(x1,...,xn) ⊕ yi = 0。
    对3×3 S盒：分量函数次数≤2，隐式方程次数自然≤2。
    对4×4 S盒：分量函数次数可达3，隐式方程次数≤3（含xi*xj*xk型）。
    返回每个方程的：次数、系数列表、可读字符串。
    """
    components = sbox_to_component_truth_tables(sbox, n, n)
    equations = []

    x_vars = [f"x{i + 1}" for i in range(n)]
    y_vars = [f"y{i + 1}" for i in range(n)]

    for i in range(n):
        anf = mobius_transform(components[i], n)
        deg = compute_algebraic_degree(anf, n)

        # 构造隐式方程字符串：fi(x) ⊕ yi = 0
        terms = []
        size = 1 << n
        for k in range(size):
            if anf[k] == 1:
                if k == 0:
                    terms.append("1")
                else:
                    factors = []
                    for bit in range(n):
                        if (k >> bit) & 1:
                            factors.append(x_vars[bit])
                    terms.append("".join(factors))
        terms.append(y_vars[i])  # 移项加上 yi

        eq_str = " ⊕ ".join(terms) + " = 0"

        # 判断该方程是否为线性（次数为1）
        is_linear = (deg <= 1)

        # 统计二次项类型（对3×3有意义）
        quadratic_terms = []
        for k in range(size):
            if anf[k] == 1 and hamming_weight(k) == 2:
                factors = []
                for bit in range(n):
                    if (k >> bit) & 1:
                        factors.append(x_vars[bit])
                quadratic_terms.append("".join(factors))

        equations.append({
            'index': i,
            'anf': anf,
            'degree': deg,
            'is_linear': is_linear,
            'equation_str': eq_str,
            'quadratic_terms': quadratic_terms,
        })

    num_linear = sum(1 for eq in equations if eq['is_linear'])

    return {
        'equations': equations,
        'num_linear_equations': num_linear,
        'implicit_class': _classify_implicit(num_linear),
    }


def _classify_implicit(num_linear):
    """按含线性方程数分类（用于3×3 S盒）"""
    if num_linear == 0:
        return 'A'
    elif num_linear == 3:
        return 'C'
    else:
        return 'B'


def print_implicit_representation(sbox, n, name="S-box", file=None):
    """格式化输出S盒的隐式代数表示"""
    result = compute_implicit_representation(sbox, n)
    lines = []
    lines.append(f"\n{'=' * 65}")
    lines.append(f"隐式代数表示：{name}")
    lines.append(f"{'=' * 65}")
    for eq in result['equations']:
        deg_label = f"（次数{eq['degree']}，{'线性' if eq['is_linear'] else '非线性'}）"
        lines.append(f"  f{eq['index'] + 1}: {eq['equation_str']}  {deg_label}")
    lines.append(f"\n  含线性方程数 k = {result['num_linear_equations']}")
    if n == 3:
        lines.append(f"  隐式分类：{result['implicit_class']}类")
    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n")
    return result