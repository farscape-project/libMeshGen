import symfem, sympy, sys

def get_basis(geom, order, dx, dy, dz):

    elem = symfem.create_element(geom, "Nedelec", order)
    basis = elem.get_basis_functions()

    if geom == "triangle":
        basis[order : 2*order] = basis[2*order - 1 : order - 1 : -1]
        basis[order : 2*order] = [-f for f in basis[order : 2*order]]
        basis = basis[2*order : 3*order] + basis[: 2*order] + basis[3*order : ]
    elif geom == "quadrilateral":
        basis[order : 2*order] = basis[2*order - 1 : order - 1 : -1]
        basis[order : 2*order] = [-f for f in basis[order : 2*order]]
        basis[3*order : 4*order] = basis[4*order - 1 : 3*order - 1 : -1]
        basis[3*order : 4*order] = [-f for f in basis[3*order : 4*order]]
        basis = basis[ : order] + basis[2*order : 4*order] + basis[order : 2*order] + basis[4*order : ]

    x, y, z = sympy.symbols('x y z')
    if geom == "quadrilateral":
        basis = [f.subs((x, y, z), ((1 + x) * .5, (1 + y) * .5, (1 + z) * .5)) for f in basis]
        basis = [f * sympy.sympify(.5) for f in basis]

    basis = [f.diff((x, dx)) for f in basis]
    basis = [f.diff((y, dy)) for f in basis]
    basis = [f.diff((z, dz)) for f in basis]

    basis = [f.with_floats() for f in basis]

    xi, eta, zeta = sympy.symbols('xi eta zeta')
    basis = [f.subs((x, y, z), (xi, eta, zeta)) for f in basis]

    for o in range(2, order + 1):
        basis = [f.subs(xi**o, sympy.UnevaluatedExpr(sympy.sympify(('xi*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs(eta**o, sympy.UnevaluatedExpr(sympy.sympify(('eta*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs(zeta**o, sympy.UnevaluatedExpr(sympy.sympify(('zeta*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]
        basis = [f.subs((xi + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(xi + 1)*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs((eta + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(eta + 1)*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs((zeta + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(zeta + 1)*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]

    return basis

dim = int(sys.argv[1])
order = int(sys.argv[2])
derivatives = int(sys.argv[3])

print("{\n"
      "  switch (elem->type())\n"
      "    {")

for geom in ["quadrilateral", "triangle"]:

    if geom == "triangle":
        print("    case TRI6:\n"
              "    case TRI7:\n"
              "      {")
    elif geom == "quadrilateral":
        print("    case QUAD8:\n"
              "    case QUAD9:\n"
              "      {")

    elem = symfem.create_reference(geom)

    if derivatives:
        print("        switch (j)\n"
              "          {")

    for d in range(derivatives + 1):
        dx = derivatives - d
        dy = d

        basis = get_basis(geom, order, dx, dy, 0)
        spaces = 6 * " " if derivatives else ""

        if derivatives:
            print(f"          case {d}:\n"
                   "            {")

        print(spaces + "        switch(ii)\n" +
              spaces + "          {")

        for f in range(len(elem.edges) * order):
            print(spaces + f"          case {f}:\n" +
                  spaces + f"            return sign * RealGradient{basis[f]};")

        for f in range(len(elem.edges) * order, len(basis)):
            print(spaces + f"          case {f}:\n" +
                  spaces + f"            return RealGradient{basis[f]};")

        print(spaces + "          default:\n" +
              spaces + "            libmesh_error_msg(\"Invalid i = \" << i);\n" +
              spaces + "          }")

        if derivatives:
            print("            }")

    if derivatives:
        print("          default:\n"
              "            libmesh_error_msg(\"Invalid j = \" << j);\n"
              "          }")

    print("      }")

print( "    default:\n"
      f"      libmesh_error_msg(\"ERROR: Unsupported {dim}D element type!: \" << Utility::enum_to_string(elem->type()));\n"
       "    }\n"
       "}")
