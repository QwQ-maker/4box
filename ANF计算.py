def mobius_transform(sbox):
    """
    输入: sbox，长度为16的列表，sbox[x] 为十六进制输出值
    输出: 4个分量布尔函数的ANF系数列表
    """
    n = 4

    # 提取各分量真值表
    # y0为最高位，y3为最低位
    components = []
    for comp in range(n):
        bit_pos = n - 1 - comp  # y0对应bit3，y3对应bit0
        f = [(sbox[x] >> bit_pos) & 1 for x in range(16)]
        components.append(f)

    # 对每个分量做Möbius变换
    anf_coeffs = []
    for f in components:
        a = f[:]
        for i in range(n):
            for x in range(16):
                if x & (1 << i):
                    a[x] ^= a[x ^ (1 << i)]
        anf_coeffs.append(a)

    return anf_coeffs


def format_anf(coeffs, n=4):
    """将系数数组格式化为ANF表达式字符串"""
    varnames = [f"x{i}" for i in range(n)]
    terms = []
    for k in range(2**n):
        if coeffs[k] == 1:
            if k == 0:
                terms.append("1")
            else:
                term = "".join(varnames[j] for j in range(n) if k & (1 << (n - 1 - j)))
                terms.append(term)
    return " ⊕ ".join(terms) if terms else "0"


def verify(sbox, anf_coeffs, x_val, n=4):
    """验证：用ANF计算x_val处的输出，与真值表对比"""
    x_bits = [(x_val >> (n - 1 - i)) & 1 for i in range(n)]
    result = []
    for comp_idx, coeffs in enumerate(anf_coeffs):
        val = 0
        for k in range(2**n):
            if coeffs[k] == 1:
                term = 1
                for j in range(n):
                    if k & (1 << (n - 1 - j)):
                        term &= x_bits[j]
                val ^= term
        result.append(val)

    computed = int("".join(map(str, result)), 2)
    expected = sbox[x_val]
    match = "✓" if computed == expected else "✗"
    print(f"  x={x_val:#04x}: ANF计算={computed:#04x}, 真值表={expected:#04x}  {match}")


# ─── S盒定义 ───────────────────────────────────────────────
sboxes = {
    "PRESENT": [0xC,0x5,0x6,0xB,0x9,0x0,0xA,0xD,0x3,0xE,0xF,0x8,0x4,0x7,0x1,0x2],
    "GIFT":    [0x1,0xA,0x4,0xC,0x6,0xF,0x3,0x9,0x2,0xD,0xB,0x7,0x5,0x0,0x8,0xE],
    "Midori":  [0xC,0xA,0xD,0x3,0xE,0xB,0xF,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6],
}

# ─── 主程序 ────────────────────────────────────────────────
for name, sbox in sboxes.items():
    print(f"\n{'='*50}")
    print(f"  {name} S盒")
    print(f"{'='*50}")

    anf = mobius_transform(sbox)

    for i, coeffs in enumerate(anf):
        expr = format_anf(coeffs)
        print(f"  y{i} = {expr}")

    print(f"\n  验证（全部16个输入）:")
    all_correct = True
    for x in range(16):
        x_bits = [(x >> (3 - i)) & 1 for i in range(4)]
        result = []
        for coeffs in anf:
            val = 0
            for k in range(16):
                if coeffs[k] == 1:
                    term = 1
                    for j in range(4):
                        if k & (1 << (3 - j)):
                            term &= x_bits[j]
                    val ^= term
            result.append(val)
        computed = int("".join(map(str, result)), 2)
        expected = sbox[x]
        if computed != expected:
            print(f"  x={x:2d}: 计算={computed:#04x}, 真值表={expected:#04x}  ✗")
            all_correct = False
    if all_correct:
        print("  所有输入验证通过 ✓")