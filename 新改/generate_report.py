"""
最终报告生成脚本：
- 16个4×4最优S盒等价类代表元的完整ANF和隐式代数表示（论文第3.6.4节用）
- 3×3双射S盒全部指纹类代表元的隐式代数表示（论文第3.6.3/5.5节用）
"""

from itertools import permutations

# ========== 基础函数 ==========

def truth_table_to_anf(tt):
    coeffs = list(tt)
    n = (len(tt)).bit_length() - 1
    for i in range(n):
        for j in range(len(tt)):
            if j & (1 << i):
                coeffs[j] ^= coeffs[j ^ (1 << i)]
    return coeffs

def anf_to_str(coeffs, n, var_names=None):
    if var_names is None:
        var_names = [f'x{i}' for i in range(n)]
    terms = []
    for u in range(1 << n):
        if coeffs[u] == 1:
            if u == 0:
                terms.append('1')
            else:
                bits = [i for i in range(n) if u & (1 << i)]
                term = ''.join(var_names[b] for b in bits)
                terms.append(term)
    return ' ⊕ '.join(terms) if terms else '0'

def anf_to_implicit_str(anf_coeffs, n, out_var):
    var_names = [f'x{i}' for i in range(n)]
    terms = []
    for u in range(1 << n):
        if anf_coeffs[u] == 1:
            if u == 0:
                terms.append('1')
            else:
                bits = [i for i in range(n) if u & (1 << i)]
                terms.append(''.join(var_names[b] for b in bits))
    terms.append(out_var)
    return ' ⊕ '.join(terms) + ' = 0'

def sbox_to_anfs(sbox, n):
    result = []
    for bit in range(n - 1, -1, -1):
        tt = [(sbox[x] >> bit) & 1 for x in range(1 << n)]
        result.append(truth_table_to_anf(tt))
    return result

def get_degree(coeffs, n):
    return max((bin(u).count('1') for u in range(1 << n) if coeffs[u] == 1), default=0)

def compute_nl(sbox, n):
    size = 1 << n
    min_nl = size
    for v in range(1, size):
        tt = [bin(v & sbox[x]).count('1') % 2 for x in range(size)]
        w = [(-1)**b for b in tt]
        step = 1
        while step < size:
            for i in range(0, size, step * 2):
                for j in range(step):
                    a, b = w[i+j], w[i+j+step]
                    w[i+j], w[i+j+step] = a+b, a-b
            step *= 2
        nl = (size // 2) - max(abs(x) for x in w) // 2
        min_nl = min(min_nl, nl)
    return min_nl

def compute_delta(sbox, n):
    size = 1 << n
    return max(
        sum(1 for x in range(size) if sbox[x] ^ sbox[x ^ dx] == dy)
        for dx in range(1, size) for dy in range(size)
    )

# ========== 16个最优4×4 S盒等价类代表元 ==========

G_reps = {
    'G0':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 11, 12,  9,  3, 14, 10,  5],
    'G1':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 11, 14,  3,  5,  9, 10, 12],
    'G2':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 11, 14,  3, 10, 12,  5,  9],
    'G3':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 12,  5,  3, 10, 14, 11,  9],
    'G4':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 12,  9, 11, 10, 14,  5,  3],
    'G5':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 12, 11,  9, 10, 14,  3,  5],
    'G6':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 12, 11,  9, 10, 14,  5,  3],
    'G7':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 12, 14, 11, 10,  9,  3,  5],
    'G8':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 14,  9,  5, 10, 11,  3, 12],
    'G9':  [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 11,  3,  5,  9, 10, 12],
    'G10': [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 11,  5, 10,  9,  3, 12],
    'G11': [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 11, 10,  5,  9, 12,  3],
    'G12': [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 11, 10,  9,  3, 12,  5],
    'G13': [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 12,  9,  5, 11, 10,  3],
    'G14': [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 12,  9, 10, 11,  5,  3],
    'G15': [0, 1, 2, 13, 4, 7, 15, 6, 8, 14, 12, 11,  3,  9,  5, 10],
}

lines = []
lines.append("=" * 80)
lines.append("【第一部分】4×4 最优S盒16个等价类代表元：显式ANF与隐式代数表示")
lines.append("数据来源：Leander & Poschmann (WAIFI 2007), Table 1")
lines.append("=" * 80)
lines.append("")
lines.append("说明：输入变量记为 x0(最高位)~x3(最低位)，输出变量 y0~y3。")
lines.append("      隐式方程由 y_i = f_i(x) 移项得 f_i(x) ⊕ y_i = 0。")
lines.append("      由于4×4 S盒分量函数最高代数次数为3，隐式方程次数≤3。")
lines.append("")

for name, sbox in G_reps.items():
    n = 4
    anf_list = sbox_to_anfs(sbox, n)
    nl = compute_nl(sbox, n)
    delta = compute_delta(sbox, n)
    degs = [get_degree(c, n) for c in anf_list]
    hex_str = ''.join(hex(v)[2:].upper() for v in sbox)
    
    lines.append(f"▶ {name}  真值表(十六进制): {hex_str}")
    lines.append(f"  密码学指标: NL={nl}, δ={delta}, deg={max(degs)},  各分量次数={degs}")
    lines.append("  显式ANF表示:")
    for i, (c, deg) in enumerate(zip(anf_list, degs)):
        anf_s = anf_to_str(c, n, [f'x{j}' for j in range(n)])
        lines.append(f"    y{i} = {anf_s}   [次数{deg}]")
    lines.append("  隐式代数表示（方程组，次数≤3）:")
    for i, c in enumerate(anf_list):
        eq = anf_to_implicit_str(c, n, f'y{i}')
        lines.append(f"    f{i}:  {eq}")
    lines.append("")

# ========== 3×3 S盒全部指纹类代表元 ==========

lines.append("")
lines.append("=" * 80)
lines.append("【第二部分】3×3 双射S盒：全部仿射指纹类代表元的隐式代数表示")
lines.append("=" * 80)
lines.append("")
lines.append("说明：3×3双射S盒共 8!=40320 个，按密码学指纹分类得到19个指纹类。")
lines.append("      输入变量 x0(最高位)~x2(最低位)，输出变量 y0~y2。")
lines.append("      由于3×3 S盒分量函数最高次数为2，隐式方程次数≤2（天然满足二次约束）。")
lines.append("")

def compute_sbox_fingerprint(sbox, n):
    nl = compute_nl(sbox, n)
    delta = compute_delta(sbox, n)
    anf_list = sbox_to_anfs(sbox, n)
    degs = tuple(sorted(get_degree(c, n) for c in anf_list))
    size = 1 << n
    lat_spec = sorted(
        sum(1 for x in range(size)
            if (bin(a & x).count('1') % 2) == (bin(b & sbox[x]).count('1') % 2)) - size // 2
        for a in range(size) for b in range(1, size)
    )
    return (nl, delta, degs, tuple(lat_spec))

print("正在枚举3×3双射S盒指纹类...")
fingerprint_map = {}
for sbox in permutations(range(8)):
    sb = list(sbox)
    fp = compute_sbox_fingerprint(sb, 3)
    if fp not in fingerprint_map:
        fingerprint_map[fp] = sb

# 统计各类S盒数量
from itertools import permutations as iperms
fp_counts = {}
for sbox in iperms(range(8)):
    fp = compute_sbox_fingerprint(list(sbox), 3)
    fp_counts[fp] = fp_counts.get(fp, 0) + 1

reps = sorted(fingerprint_map.items(), key=lambda x: (x[0][1], -x[0][0]))

for idx, (fp, rep) in enumerate(reps):
    nl_v, delta_v, degs_sorted, _ = fp
    n = 3
    anf_list = sbox_to_anfs(rep, n)
    degs = [get_degree(c, n) for c in anf_list]
    linear_count = sum(1 for d in degs if d <= 1)
    cls = 'A类(全二次)' if linear_count == 0 else ('C类(全线性)' if linear_count == 3 else 'B类(含线性)')
    count = fp_counts[fp]
    hex_str = ''.join(hex(v)[2:].upper() for v in rep)
    
    lines.append(f"类 #{idx+1}  真值表: {rep}  (十六进制: {hex_str})")
    lines.append(f"  密码学指标: NL={nl_v}, δ={delta_v}, deg={max(degs)},  各分量次数={degs}")
    lines.append(f"  本类S盒总数: {count}  隐式结构类型: {cls}（线性方程数 k={linear_count}）")
    lines.append("  显式ANF表示:")
    for i, c in enumerate(anf_list):
        s = anf_to_str(c, n, [f'x{j}' for j in range(n)])
        lines.append(f"    y{i} = {s}   [次数{degs[i]}]")
    lines.append("  隐式代数表示（方程组，次数≤2）:")
    for i, c in enumerate(anf_list):
        eq = anf_to_implicit_str(c, n, f'y{i}')
        lines.append(f"    f{i}:  {eq}")
    lines.append("")

# 统计汇总
lines.append("=" * 80)
lines.append("【汇总统计】3×3 双射S盒分类汇总")
lines.append("=" * 80)
lines.append("")
from collections import defaultdict
stat = defaultdict(int)
for fp, cnt in fp_counts.items():
    nl_v, delta_v, degs_s, _ = fp
    stat[(nl_v, delta_v)] += cnt
lines.append(f"{'(NL, δ)':<15}{'S盒数量':<12}{'占比':<10}密码学评价")
lines.append("-" * 60)
total = 40320
for (nl_v, delta_v), cnt in sorted(stat.items(), key=lambda x: (-x[0][0], x[0][1])):
    pct = cnt / total * 100
    eval_s = '最优' if nl_v == 2 and delta_v == 2 else ('较差' if delta_v == 8 else '一般')
    lines.append(f"({nl_v}, {delta_v})         {cnt:<12}{pct:.2f}%     {eval_s}")
lines.append("-" * 60)
lines.append(f"{'合计':<15}{total:<12}{100:.2f}%")
lines.append("")
lines.append("关键结论：NL=2与δ=2对于3×3双射S盒是等价条件（两者同时成立或同时不成立）。")
lines.append("最优S盒（NL=2, δ=2）共10752个，占26.67%，均属A类（全二次隐式方程）。")

report = '\n'.join(lines)
with open('sbox_report.txt', 'w', encoding='utf-8') as f:
    f.write(report)
print("报告已生成：sbox_report.txt")
print(f"总行数：{len(lines)}")
