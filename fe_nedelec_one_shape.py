import inflect, itertools, re, symfem, sympy, sys

def get_basis(geom, order, dx, dy, dz):

    elem = symfem.create_element(geom, "Nedelec", order)
    basis = elem.get_basis_functions()

    pe = len(elem.entity_dofs(1, 0))
    fo = pe * len(elem.reference.sub_entities(1))
    pf = len(elem.entity_dofs(2, 0))
    vo = fo + pf * len(elem.reference.sub_entities(2))

    se  = lambda edge : slice(pe * edge, pe * (edge + 1))
    sre = lambda edge : slice(pe * (edge + 1) - 1, pe * edge - 1 if edge > 0 else None, -1)
    sf  = lambda face : slice(fo + pf * face, fo + pf * (face + 1))
    srf = lambda face : slice(fo + pf * (face + 1) - 1, fo + pf * face - 1, -1)

    for e in range(len(elem.reference.edges)):
        edge_dofs = elem.dof_plot_positions()[se(e)]
        edge_funs = elem.get_basis_functions()[se(e)]
        basis[se(e)] = [f for _, f in sorted(zip(edge_dofs, edge_funs), key = lambda x: x[0])]

    if geom == "triangle":
        basis[se(0)] = basis[sre(0)]
        basis[se(1)] = [-f for f in basis[sre(1)]]
        basis = [b for e in (2, 0, 1) for b in basis[se(e)]] + basis[fo : ]
    elif geom == "quadrilateral":
        basis[se(1)] = [-f for f in basis[sre(1)]]
        basis[se(3)] = [-f for f in basis[sre(3)]]
        basis = [b for e in (0, 2, 3, 1) for b in basis[se(e)]] + basis[fo : ]
    elif geom == "tetrahedron":
        basis[se(0)] = basis[sre(0)]
        basis[se(1)] = basis[sre(1)]
        basis[se(2)] = basis[sre(2)]
        basis = [b for e in (5, 2, 4, 3, 1, 0) for b in basis[se(e)]] +  [b for f in (3, 2, 0, 1) for b in basis[sf(f)]] + basis[vo : ]
    elif geom == "hexahedron":
        basis[se(5)] = [-f for f in basis[sre(5)]]
        basis[se(11)] = [-f for f in basis[sre(11)]]
        basis = [b for e in (0, 3, 5, 1, 2, 4, 7, 6, 8, 10, 11, 9) for b in basis[se(e)]] +  [b for f in (0, 1, 3, 4, 2, 5) for b in basis[sf(f)]] + basis[vo : ]

    x, y, z = sympy.symbols('x y z')
    if geom == "quadrilateral" or geom == "hexahedron":
        basis = [f.subs((x, y, z), ((1 + x) / 2, (1 + y) / 2, (1 + z) / 2)) for f in basis]
        basis = [f / sympy.sympify(2) for f in basis]

    for v, d in ((x, dx), (y, dy), (z, dz)):
        basis = [f.diff((v, d)) for f in basis]

    xi, eta, zeta = sympy.symbols('xi eta zeta')
    basis = [f.subs((x, y, z), (xi, eta, zeta)) for f in basis]

    for o in range(order, 1, -1):
        for vsym, vstr in ((xi, 'xi'), (eta, 'eta'), (zeta, 'zeta')):
            for esym, estr in ((vsym**o, (vstr+'*')*o), ((vsym + 1)**o, ('('+vstr+' + 1)*')*o), (((vsym + 1)/2)**o, ('('+vstr+' + 1)/2*')*o)):
                basis = [f.subs(esym, sympy.UnevaluatedExpr(sympy.sympify((estr)[:-1], locals={vstr: vsym}, evaluate = False))) for f in basis]

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
