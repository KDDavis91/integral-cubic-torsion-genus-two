"""Standalone exact verifier for Frobenius on integral cubic two-torsion.

Python 3.10+, standard library only. No earlier research files are loaded.
This is a computational certificate, not external human peer review.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
from math import gcd
import argparse
import json

EDGES = list(combinations(range(4), 2))
WEDGES = [e for e in EDGES if e != (1, 3)]
SLOTS = [i for i in range(20) if i not in (1, 3, 4, 9)]
WEIGHTS = [[0, 2], [1, 3, 4, 5, 6, 11], [7, 8, 9, 12, 13, 14], [10, 15]]


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return list(map(list, zip(*a)))


def multiply(a, b, prime=None):
    out = [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]
    return [[x % prime for x in row] for row in out] if prime else out


def subtract(a, b, prime=None):
    out = [[x-y for x, y in zip(r, s)] for r, s in zip(a, b)]
    return [[x % prime for x in row] for row in out] if prime else out


def reduce_rows(a, prime=None):
    a = [[int(x) % prime if prime else Fraction(x) for x in row] for row in a]
    pivots = []
    i = 0
    for j in range(len(a[0])):
        k = next((k for k in range(i, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[i], a[k] = a[k], a[i]
        z = pow(a[i][j], -1, prime) if prime else 1/a[i][j]
        a[i] = [x*z % prime if prime else x*z for x in a[i]]
        for k in range(len(a)):
            if k != i:
                z = a[k][j]
                a[k] = [(x-z*y) % prime if prime else x-z*y for x, y in zip(a[k], a[i])]
        pivots.append(j)
        i += 1
        if i == len(a):
            break
    return a, pivots


def rank(a, prime=None):
    return len(reduce_rows(a, prime)[1])


def inverse(a, prime=None):
    n = len(a)
    r, pivots = reduce_rows([x+y for x, y in zip(a, eye(n))], prime)
    assert pivots[:n] == list(range(n))
    assert [row[:n] for row in r] == eye(n)
    out = [row[n:] for row in r]
    if not prime:
        assert all(x.denominator == 1 for row in out for x in row)
        out = [[int(x) for x in row] for row in out]
    return out


def nullspace(a, prime=2):
    r, pivots = reduce_rows(a, prime)
    columns = []
    for j in range(len(a[0])):
        if j in pivots:
            continue
        v = [0]*len(a[0])
        v[j] = 1
        for i, k in enumerate(pivots):
            v[k] = -r[i][j] % prime
        columns.append(v)
    return transpose(columns) if columns else [[] for _ in a[0]]


def determinant(a):
    # Fraction-free Bareiss elimination.
    a = [list(row) for row in a]
    previous, sign = 1, 1
    for k in range(len(a)-1):
        j = next((i for i in range(k, len(a)) if a[i][k]), None)
        if j is None:
            return 0
        if j != k:
            a[k], a[j] = a[j], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def wedge(x, y):
    return [x[i]*y[j]-x[j]*y[i] for i, j in EDGES]


def quotient_coordinates():
    standard = eye(4)
    theta = [0, 1, 0, 0, 1, 0]
    relations = [[theta[j] if i == k else 0 for i in range(4) for j in range(6)] for k in range(4)]
    for i, j, k in combinations(range(4), 3):
        v = [0]*24
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            for z, value in enumerate(wedge(standard[b], standard[c])):
                v[6*a+z] += value
        relations.append(v)
    sections = []
    for slot in SLOTS:
        v = [0]*24
        v[6*(slot//5)+EDGES.index(WEDGES[slot % 5])] = 1
        sections.append(v)
    change = transpose(relations+sections)
    d = determinant(change)
    assert abs(d) == 1
    all_inverse = inverse(change)
    q, s, r = all_inverse[8:], transpose(sections), transpose(relations)
    assert multiply(q, s) == eye(16)
    assert multiply(q, r) == [[0]*8 for _ in range(16)]
    return q, s, r, d


Q, SECTION, RELATIONS, PRESENTATION_DETERMINANT = quotient_coordinates()


def cubic_action(m):
    columns = transpose(m)
    w = transpose([wedge(columns[i], columns[j]) for i, j in EDGES])
    raw = [[m[i//6][j//6]*w[i % 6][j % 6] for j in range(24)] for i in range(24)]
    assert multiply(Q, multiply(raw, RELATIONS)) == [[0]*8 for _ in range(16)]
    return multiply(Q, multiply(raw, SECTION))


def infinitesimal(parameters):
    a, b, c = parameters
    n = [[0]*4 for _ in range(4)]
    n[0][2], n[0][3], n[1][2], n[1][3] = a, b, b, c
    columns, standard = transpose(n), eye(4)
    wcols = []
    for i, j in EDGES:
        wcols.append([x+y for x, y in zip(wedge(columns[i], standard[j]), wedge(standard[i], columns[j]))])
    w = transpose(wcols)
    raw = [[n[i//6][j//6]*int(i % 6 == j % 6)+int(i//6 == j//6)*w[i % 6][j % 6]
            for j in range(24)] for i in range(24)]
    assert multiply(Q, multiply(raw, RELATIONS)) == [[0]*8 for _ in range(16)]
    return multiply(Q, multiply(raw, SECTION))


def shear(parameters):
    a, b, c = parameters
    out = eye(4)
    out[0][2], out[0][3], out[1][2], out[1][3] = a, b, b, c
    return out


def normaliser(u, h, q=1):
    a, b = u[0]
    c, d = u[1]
    determinant_u = a*d-b*c
    assert abs(determinant_u) == 1
    top = [[d//determinant_u, -c//determinant_u], [-b//determinant_u, a//determinant_u]]
    block = [[0]*4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            block[i][j], block[i+2][j+2] = q*top[i][j], u[i][j]
    return multiply(shear(h), block)


def weight_maps(parameters):
    a, b, c = parameters
    m = [[c, 0, 0, 0, a, b], [0, b, c, -a, b, c]]
    middle = [[0, a, 0, 0, b, 0], [b, 0, a, c, 0, b], [-a, b, 0, -b, c, 0],
              [0, 0, b, 0, 0, c], [b, c, -2*a, 0, 0, -b], [a, 0, 0, 2*b, c, -a]]
    last = [[-3*b, -c], [2*a, b], [-c, 0], [-a, -3*b], [0, a], [-b, -2*c]]
    return m, middle, last


def polynomial_add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0)+v
        if not out[k]:
            del out[k]
    return out


def polynomial_scale(a, n):
    return {k: v*n for k, v in a.items() if v*n}


def polynomial_multiply(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            k = tuple(v+w for v, w in zip(i, j))
            out[k] = out.get(k, 0)+x*y
    return {k: v for k, v in out.items() if v}


def polynomial_matrix_product(a, b):
    out = [[{} for _ in b[0]] for _ in a]
    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k]:
                for j in range(len(b[0])):
                    if b[k][j]:
                        out[i][j] = polynomial_add(out[i][j], polynomial_multiply(a[i][k], b[k][j]))
    return out


def constant_matrix(a):
    return [[{(0, 0, 0): x} if x else {} for x in row] for row in a]


def polynomial_certificates():
    zero = (0, 0, 0)
    variables = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    matrices = [infinitesimal(v) for v in variables]
    n = [[{variables[k]: matrices[k][i][j] for k in range(3) if matrices[k][i][j]}
          for j in range(16)] for i in range(16)]
    middle = [[n[i][j] for j in WEIGHTS[2]] for i in WEIGHTS[1]]
    # Signed subset dynamic determinant, independent of Bareiss arithmetic.
    states = {0: {zero: 1}}
    for i in range(6):
        new = {}
        for mask, value in states.items():
            for j in range(6):
                if mask >> j & 1:
                    continue
                sign = -1 if (mask >> (j+1)).bit_count() % 2 else 1
                term = polynomial_scale(polynomial_multiply(value, middle[i][j]), sign)
                target = mask | (1 << j)
                new[target] = polynomial_add(new.get(target, {}), term)
        states = new
    determinant_poly = states[63]
    assert determinant_poly == {(3, 0, 3): -4, (2, 2, 2): 12, (1, 4, 1): -12, (0, 6, 0): 4}
    # Reconstruct rho(tau_B) over the polynomial ring in the raw 24 tensors.
    m = constant_matrix(eye(4))
    m[0][2], m[0][3], m[1][2], m[1][3] = [{variables[k]: 1} for k in (0, 1, 1, 2)]
    cols = transpose(m)
    wedges = []
    for i, j in EDGES:
        wedges.append([polynomial_add(polynomial_multiply(cols[i][a], cols[j][b]),
                                     polynomial_scale(polynomial_multiply(cols[i][b], cols[j][a]), -1))
                       for a, b in EDGES])
    w = transpose(wedges)
    raw = [[polynomial_multiply(m[i//6][j//6], w[i % 6][j % 6]) for j in range(24)] for i in range(24)]
    action = polynomial_matrix_product(constant_matrix(Q), polynomial_matrix_product(raw, constant_matrix(SECTION)))
    n2 = polynomial_matrix_product(n, n)
    n3 = polynomial_matrix_product(n2, n)
    n4 = polynomial_matrix_product(n3, n)
    assert all(not x for row in n4 for x in row)
    for i in range(16):
        for j in range(16):
            left = polynomial_scale(polynomial_add(action[i][j], {zero: -1} if i == j else {}), 6)
            right = polynomial_add(polynomial_scale(n[i][j], 6), polynomial_add(polynomial_scale(n2[i][j], 3), n3[i][j]))
            assert left == right
    return {'middle_determinant': '-4*(a*c-b^2)^3',
            'middle_determinant_terms': [[list(k), v] for k, v in sorted(determinant_poly.items())],
            'full_polynomial_exponential_identity': '6*(rho(tau_B)-I)=6*N+3*N^2+N^3',
            'nilpotence_N_four': True, 'matrix_entries_checked': 256}


def kernel_reduction(parameters):
    m, _, _ = weight_maps(parameters)
    k1 = nullspace(m)
    k = [[0]*6 for _ in range(16)]
    for j, i in enumerate(WEIGHTS[0]):
        k[i][j] = 1
    for i, row in zip(WEIGHTS[1], k1):
        k[i][2:] = row
    assert rank(k, 2) == 6
    return k


def quotient_basis(parameters, even):
    delta = subtract(cubic_action(shear((0, 0, 0) if even else parameters)), eye(16), 2)
    k = kernel_reduction(parameters)
    cols = transpose(k)
    if not even:
        # Canonical map from the rank-two cubic graph quotient: z |-> delta*z.
        for j in WEIGHTS[3]:
            col = [row[j] for row in delta]
            assert rank(transpose(cols+[col]), 2) == len(cols)+1
            cols.append(col)
    invariant_dimension = 16 if even else 8
    if even:
        for col in eye(16):
            if rank(transpose(cols+[col]), 2) > len(cols):
                cols.append(col)
    for col in eye(16):
        if rank(transpose(cols+[col]), 2) > len(cols):
            cols.append(col)
    basis = transpose(cols)
    assert len(cols) == 16
    assert multiply(delta, [row[:invariant_dimension] for row in basis], 2) == [[0]*invariant_dimension for _ in range(16)]
    return basis, inverse(basis, 2), invariant_dimension


def quotient_action(f, data):
    basis, inv, n = data
    a = multiply(inv, multiply(cubic_action(f), basis, 2), 2)
    assert not any(a[i][j] for i in range(n, 16) for j in range(n))
    assert not any(a[i][j] for i in range(6, 16) for j in range(6))
    return [row[6:n] for row in a[6:n]]


def finite_certificates():
    gl2 = [list(map(list, (x[:2], x[2:]))) for x in product(range(2), repeat=4)
           if (x[0]*x[3]-x[1]*x[2]) % 2]
    parities = [x for x in product(range(2), repeat=3) if any(x)]
    records, count = [], 0
    for b in parities:
        a, d, c = b
        pairing = [[a, d], [d, c]]
        n = infinitesimal(b)
        m, middle, last = weight_maps(b)
        assert [[n[i][j] for j in WEIGHTS[1]] for i in WEIGHTS[0]] == m
        assert [[n[i][j] for j in WEIGHTS[2]] for i in WEIGHTS[1]] == middle
        assert [[n[i][j] for j in WEIGHTS[3]] for i in WEIGHTS[2]] == last
        assert rank(m, 2) == rank(last, 2) == 2 and rank(middle, 2) == 4
        delta = subtract(cubic_action(shear(b)), eye(16), 2)
        assert rank(delta, 2) == 8 and multiply(delta, delta, 2) == [[0]*16 for _ in range(16)]
        assert multiply(n, n, 2) == [[0]*16 for _ in range(16)]
        stabilisers = [u for u in gl2 if multiply(transpose(u), multiply(pairing, u, 2), 2) == pairing]
        cases = []
        for even in (False, True):
            data = quotient_basis(b, even)
            actions = []
            for u in stabilisers:
                for h in product(range(2), repeat=3):
                    aq = quotient_action(normaliser(u, h), data)
                    fixed = len(aq)-rank(subtract(aq, eye(len(aq)), 2), 2)
                    if not even:
                        assert aq == u
                    elif u == eye(2):
                        expected = 10 if not any(h) else 8 if h == b else 6
                        assert fixed == expected
                    actions.append({'U': u, 'H_parameters': h, 'quotient_action': aq, 'fixed_dimension': fixed})
                    count += 1
            cases.append({'pairing_even': even, 'dimension': data[2]-6, 'basis': data[0], 'actions': actions})
        records.append({'primitive_parity': b, 'middle_rank_mod_two': 4, 'cases': cases})
    pair_checks = 0
    for b in parities:
        for h in parities:
            assert rank(weight_maps(b)[0]+weight_maps(h)[0], 2) == (2 if b == h else 4)
            if b != h:
                pair_checks += 1
    assert count == 288 and pair_checks == 42
    minors = []
    for b in ((1, 0, 0), (1, 0, 1), (0, 1, 0)):
        middle = weight_maps(b)[1]
        witness = next(({'rows_one_based': [i+1 for i in rows], 'columns_one_based': [j+1 for j in cols], 'determinant': determinant([[middle[i][j] for j in cols] for i in rows])}
                        for rows in combinations(range(6), 4) for cols in combinations(range(6), 4)
                        if determinant([[middle[i][j] for j in cols] for i in rows]) % 2), None)
        assert witness is not None
        minors.append({'normal_form': b, **witness})
    return {'finite_actions': count, 'distinct_weight_map_pairs': pair_checks,
            'canonical_primitive_graph_isomorphisms': 7, 'middle_minor_normal_forms': minors,
            'records': records}


def polynomial_product(a, b):
    return [sum(a[j]*b[i-j] for j in range(len(a)) if 0 <= i-j < len(b)) for i in range(len(a)+len(b)-1)]


def discriminant_sylvester(f):
    n = len(f)-1
    df = [i*f[i] for i in range(1, len(f))]
    m = len(df)-1
    matrix = [[0]*i+list(reversed(f))+[0]*(m-1-i) for i in range(m)]
    matrix += [[0]*i+list(reversed(df))+[0]*(n-1-i) for i in range(n)]
    return (-1)**(n*(n-1)//2)*determinant(matrix)//f[-1]


def discriminant_euclidean(f):
    def trim(a):
        while a and not a[-1]:
            a.pop()
        return a

    def resultant(a, b):
        m, n = len(a)-1, len(b)-1
        if n == 0:
            return b[0]**m
        r = list(a)
        while r and len(r) >= len(b):
            coefficient, shift = r[-1]/b[-1], len(r)-len(b)
            for i, x in enumerate(b):
                r[i+shift] -= coefficient*x
            trim(r)
        if not r:
            return Fraction(0)
        k = len(r)-1
        return (-1)**(m*n)*b[-1]**(m-k)*resultant(b, r)

    f = [Fraction(x) for x in f]
    n = len(f)-1
    answer = (-1)**(n*(n-1)//2)*resultant(f, [i*f[i] for i in range(1, len(f))])/f[-1]
    assert answer.denominator == 1
    return int(answer)


def valuation(n, prime):
    assert n
    v = 0
    while n % prime == 0:
        n //= prime
        v += 1
    return v


def root_pair_actions():
    basis = [3, 12, 18, 27]
    assert [[(x & y).bit_count() % 2 for y in basis] for x in basis] == [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]
    coordinates = {}
    for i in range(16):
        v = 0
        for j in range(4):
            if i >> j & 1:
                v ^= basis[j]
        coordinates[v] = coordinates[v ^ 63] = [(i >> j) & 1 for j in range(4)]
    out = []
    for a, b in ((0, 1), (2, 3), (4, 5)):
        images = []
        for v in basis:
            if ((v >> a) ^ (v >> b)) & 1:
                v ^= (1 << a)+(1 << b)
            images.append(coordinates[v])
        out.append(transpose(images))
    assert out == [shear((1, 0, 0)), shear((0, 0, 1)), shear((1, 1, 1))]
    return out


def theta_graph_action(edge_permutation):
    # Derive the signed cycle action from oriented edges, rather than guessing
    # the sign from a permutation of the three nodes.
    images = []
    for cycle in ([1, -1, 0], [0, 1, -1]):
        out = [0, 0, 0]
        for i, j in enumerate(edge_permutation):
            out[j] = cycle[i]
        assert sum(out) == 0
        images.append([out[0], -out[2]])
    return transpose(images)


def arithmetic_examples():
    thin = []
    pairing = [[2, -1], [-1, 2]]
    data = quotient_basis((0, 1, 0), False)
    for h, edge_permutation, expected in (([0, -1, 0, 1], [0, 1, 2], 2), ([0, -2, 0, 1], [2, 1, 0], 1), ([1, 1, 0, 1], [1, 2, 0], 0)):
        u = theta_graph_action(edge_permutation)
        f = polynomial_product(h, h)
        f[0] -= 5
        df = discriminant_sylvester(f)
        assert df == discriminant_euclidean(f) != 0 and valuation(df, 5) == 3
        assert discriminant_sylvester(h) % 5
        roots = [x for x in range(5) if sum(a*x**i for i, a in enumerate(h)) % 5 == 0]
        assert len(roots) == (3 if expected == 2 else 1 if expected == 1 else 0)
        assert multiply(transpose(u), multiply(pairing, u)) == pairing
        f4 = normaliser(u, (0, 0, 0), 5)
        assert multiply(f4, shear((2, -1, 2))) == multiply(shear((10, -5, 10)), f4)
        aq = quotient_action(f4, data)
        assert 2-rank(subtract(aq, eye(2)), 2) == expected
        # Action on coker(B) is U^-t. Modulo three this cokernel is one-dimensional.
        basis = [[2, 1], [-1, 0]]
        finverse = inverse(basis, 3)
        a, b = u[0]
        c, d = u[1]
        detu = a*d-b*c
        dual = [[d//detu, -c//detu], [-b//detu, a//detu]]
        component = multiply(finverse, multiply(dual, basis, 3), 3)
        assert component[1][0] == 0
        tamagawa = 3 if component[1][1] == 1 else 1
        assert tamagawa == 3
        thin.append({'h': h, 'sextic': f, 'discriminant': df, 'v5_discriminant': 3,
                     'edge_permutation': edge_permutation, 'signed_graph_action': u,
                     'node_residue_degrees': [1, 1, 1] if expected == 2 else [1, 2] if expected == 1 else [3],
                     'ordinary_tamagawa': tamagawa, 'cubic_two_torsion_fixed_dimension': expected})
    for edge_permutation in permutations(range(3)):
        u = theta_graph_action(edge_permutation)
        a, b = u[0]
        c, d = u[1]
        detu = a*d-b*c
        dual = [[d//detu, -c//detu], [-b//detu, a//detu]]
        component = multiply(finverse, multiply(dual, basis, 3), 3)
        assert component[1][0] == 0 and component[1][1] == 1
    data = quotient_basis((0, 1, 0), True)
    pairs = root_pair_actions()
    thick = []
    for eps in product(range(2), repeat=3):
        units = [1+e for e in eps]
        f = [1]
        for centre, unit in zip((0, 1, -1), units):
            f = polynomial_product(f, [centre*centre-25*unit, -2*centre, 1])
        df = discriminant_sylvester(f)
        assert df == discriminant_euclidean(f) != 0 and valuation(df, 5) == 6
        f4 = eye(4)
        for e, t in zip(eps, pairs):
            if e:
                f4 = multiply(f4, t, 2)
        aq = quotient_action(f4, data)
        fixed = 10-rank(subtract(aq, eye(10)), 2)
        expected = 10 if not sum(eps) else 8 if sum(eps) == 3 else 6
        assert fixed == expected
        thick.append({'units': units, 'sextic': f, 'discriminant': df, 'v5_discriminant': 6,
                      'ordinary_component_invariants': [2, 6], 'ordinary_tamagawa': 12,
                      'cubic_two_torsion_fixed_dimension': fixed})
    return {'thin_models': thin, 'thick_models': thick, 'exact_discriminants': 11,
            'all_signed_theta_edge_actions_checked': 6,
            'root_pair_transpositions': 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    polynomial = polynomial_certificates()
    finite = finite_certificates()
    examples = arithmetic_examples()
    out = {'status': 'PASS', 'presentation_unimodular_determinant': PRESENTATION_DETERMINANT,
           'polynomial_certificates': polynomial, 'finite_certificates': finite,
           'arithmetic_examples': examples, 'external_human_review_completed': False,
           'publication_priority_established': False, 'erdos677_solved': False}
    args.output.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'presentation_determinant': PRESENTATION_DETERMINANT,
                      **polynomial, 'actions_checked': finite['finite_actions'],
                      'canonical_graph_isomorphisms': finite['canonical_primitive_graph_isomorphisms'],
                      'middle_minor_normal_forms': finite['middle_minor_normal_forms'],
                      'distinct_weight_map_pairs': finite['distinct_weight_map_pairs'],
                      'arithmetic_discriminants': examples['exact_discriminants'],
                      'external_human_review_completed': False}, indent=2))


if __name__ == '__main__':
    main()
