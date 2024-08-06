import inflect, itertools, re, symfem, sympy, sys

def get_basis(geom, order, dx, dy, dz):

    elem = symfem.create_element(geom, "Nedelec", order)
    basis = elem.get_basis_functions()

    pe = order
    fo = pe * len(elem.reference.edges)
    pf = order*(order - 1) if geom == "tetrahedron" else 2*order*(order - 1) if geom == "hexahedron" else 0

    for e in range(len(elem.reference.edges)):
        edge_dofs = elem.dof_plot_positions()[e*pe : (e+1)*pe]
        edge_funs = elem.get_basis_functions()[e*pe : (e+1)*pe]
        basis[e*pe : (e+1)*pe] = [f for _, f in sorted(zip(edge_dofs, edge_funs), key = lambda x: x[0])]

    if geom == "triangle":
        basis[: pe] = basis[pe - 1 : : -1]
        basis[pe : 2*pe] = basis[2*pe - 1 : pe - 1 : -1]
        basis[pe : 2*pe] = [-f for f in basis[pe : 2*pe]]
        basis = basis[2*pe : 3*pe] + basis[: 2*pe] + basis[fo : ]
    elif geom == "quadrilateral":
        basis[pe : 2*pe] = basis[2*pe - 1 : pe - 1 : -1]
        basis[pe : 2*pe] = [-f for f in basis[pe : 2*pe]]
        basis[3*pe : 4*pe] = basis[4*pe - 1 : 3*pe - 1 : -1]
        basis[3*pe : 4*pe] = [-f for f in basis[3*pe : 4*pe]]
        basis = basis[ : pe] + basis[2*pe : 4*pe] + basis[pe : 2*pe] + basis[fo : ]
    elif geom == "tetrahedron":
        basis[ : pe] = basis[pe - 1 : : -1]
        basis[pe : 2*pe] = basis[2*pe - 1 : pe - 1 : -1]
        basis[2*pe : 3*pe] = basis[3*pe - 1 : 2*pe - 1 : -1]
        basis = basis[5*pe : 6*pe] + basis[2*pe : 3*pe] + basis[4*pe : 5*pe] + basis[3*pe : 4*pe] + basis[pe : 2*pe] + basis[ : pe] + basis[fo : ]
        basis = basis[: fo] + basis[fo + 3*pf : fo + 4*pf] + basis[fo + 2*pf : fo + 3*pf] + basis[fo : fo + 2*pf] + basis[fo + 4*pf : ]
    elif geom == "hexahedron":
        basis[5*pe : 6*pe] = basis[6*pe - 1 : 5*pe - 1 : -1]
        basis[5*pe : 6*pe] = [-f for f in basis[5*pe : 6*pe]]
        basis[11*pe : 12*pe] = basis[12*pe - 1 : 11*pe - 1 : -1]
        basis[11*pe : 12*pe] = [-f for f in basis[11*pe : 12*pe]]
        basis = basis[ : pe] + basis[3*pe : 4*pe] + basis[5*pe : 6*pe] + basis[pe : 3*pe] + basis[4*pe : 5*pe] + basis[7*pe : 8*pe] + basis[6*pe : 7*pe] + basis[8*pe : 9*pe] + basis[10*pe : 12*pe] + basis[9*pe : 10*pe] + basis[fo : ]
        basis = basis[: fo + 2*pf] + basis[fo + 3*pf : fo + 5*pf] + basis[fo + 2*pf : fo + 3*pf] + basis[fo + 5*pf : ]

    x, y, z = sympy.symbols('x y z')
    if geom == "quadrilateral" or geom == "hexahedron":
        basis = [f.subs((x, y, z), ((1 + x) / 2, (1 + y) / 2, (1 + z) / 2)) for f in basis]
        basis = [f / sympy.sympify(2) for f in basis]

    basis = [f.diff((x, dx)) for f in basis]
    basis = [f.diff((y, dy)) for f in basis]
    basis = [f.diff((z, dz)) for f in basis]

    xi, eta, zeta = sympy.symbols('xi eta zeta')
    basis = [f.subs((x, y, z), (xi, eta, zeta)) for f in basis]

    for o in range(order, 1, -1):
        basis = [f.subs(xi**o, sympy.UnevaluatedExpr(sympy.sympify(('xi*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs(eta**o, sympy.UnevaluatedExpr(sympy.sympify(('eta*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs(zeta**o, sympy.UnevaluatedExpr(sympy.sympify(('zeta*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]
        basis = [f.subs((xi + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(xi + 1)*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs((eta + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(eta + 1)*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs((zeta + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(zeta + 1)*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]
        basis = [f.subs(((xi + 1)/2)**o, sympy.UnevaluatedExpr(sympy.sympify(('(xi + 1)/2*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs(((eta + 1)/2)**o, sympy.UnevaluatedExpr(sympy.sympify(('(eta + 1)/2*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs(((zeta + 1)/2)**o, sympy.UnevaluatedExpr(sympy.sympify(('(zeta + 1)/2*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]

    p = re.compile(r'(\d+)')
    basis = [p.sub(r'\1.', str(f)) for f in basis]

    return basis

dim = int(sys.argv[1])
order = int(sys.argv[2])
derivatives = int(sys.argv[3])

p = inflect.engine()

print("case " + p.number_to_words(p.ordinal(order)).upper() + ":\n"
      "  {\n"
      "    switch (elem->type())\n"
      "      {")

for geom in ["quadrilateral", "triangle"] if dim == 2 else \
            ["hexahedron", "tetrahedron"] if dim == 3 else \
            []:

    if geom == "triangle":
        print("      case TRI6:\n"
              "      case TRI7:\n"
              "        {")
    elif geom == "quadrilateral":
        print("      case QUAD8:\n"
              "      case QUAD9:\n"
              "        {")
    elif geom == "tetrahedron":
        print("      case TET10:") if order < 2 else None
        print("      case TET14:\n"
              "        {")
    elif geom == "hexahedron":
        print("      case HEX20:") if order < 2 else None
        print("      case HEX27:\n"
              "        {")

    elem = symfem.create_reference(geom)

    if derivatives:
        print("          switch (j)\n"
              "            {")

    combs = []
    for d in itertools.combinations(range(derivatives + dim - 1), dim - 1):
        combs.append([b - a - 1 for a, b in zip((-1,) + d, d + (derivatives + dim - 1,))])
    combs = combs[::-1]

    if dim == 2:
        for d in combs: d.append(0)
    elif dim == 3 and derivatives == 2:
        combs[2], combs[3] = combs[3], combs[2]

    for d in range(len(combs)):
        dx, dy, dz = combs[d]

        basis = get_basis(geom, order, dx, dy, dz)
        spaces = 6 * " " if derivatives else ""

        if derivatives:
            print(f"              // d" + f"^{derivatives}" * (derivatives > 1) + "()/" +
                                    "dxi" * (dx > 0) + f"^{dx}" * (dx > 1) +
                                    "deta" * (dy > 0) + f"^{dy}" * (dy > 1) +
                                    "dzeta" * (dz > 0) + f"^{dz}" * (dz > 1) + "\n"
                  f"            case {d}:\n"
                   "              {")

        print(spaces + "          switch(ii)\n" +
              spaces + "            {")

        for f in range(len(elem.edges) * order):
            print(spaces + f"            case {f}:\n" +
                  spaces + f"              return sign * RealGradient{basis[f]};")

        for f in range(len(elem.edges) * order, len(basis)):
            print(spaces + f"            case {f}:\n" +
                  spaces + f"              return RealGradient{basis[f]};")

        print(spaces + "            default:\n" +
              spaces + "              libmesh_error_msg(\"Invalid i = \" << i);\n" +
              spaces + "            }")

        if derivatives:
            print(f"              }} // j = {d}\n")

    if derivatives:
        print("            default:\n"
              "              libmesh_error_msg(\"Invalid j = \" << j);\n"
              "            }")

    print("        }\n")

print( "      default:\n"
      f"        libmesh_error_msg(\"ERROR: Unsupported {dim}D element type!: \" << Utility::enum_to_string(elem->type()));\n"
       "      } // end switch (type)\n"
       "  } // end case " + p.number_to_words(p.ordinal(order)).upper())
