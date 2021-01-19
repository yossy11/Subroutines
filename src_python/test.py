import sympy
F, G, H, L, M, N = sympy.symbols("F G H L M N")
sxx, syy, szz, sxy, sxz, syz = sympy.symbols("sxx syy szz sxy sxz syz")
eq_stress = sympy.sqrt(1.5*(F*(syy-szz)**2 + G*(szz-sxx)**2 +
                            H*(sxx-syy) ** 2 + 2*L*syz**2 + 2*M*sxz**2 + 2*N*sxy**2)/(F+G+H))
dgdsxx = sympy.diff(eq_stress, sxx)
dgdsyy = sympy.diff(eq_stress, syy)
dgdszz = sympy.diff(eq_stress, szz)
dgdsxy = sympy.diff(eq_stress, sxy)
dgdsxz = sympy.diff(eq_stress, sxz)
dgdsyz = sympy.diff(eq_stress, syz)

ddgdsxxdsxx = sympy.diff(dgdsxx, sxx)
ddgdsxxdsyy = sympy.diff(dgdsxx, syy)
ddgdsxxdszz = sympy.diff(dgdsxx, szz)
ddgdsxxdsxy = sympy.diff(dgdsxx, sxy)
ddgdsxxdsxz = sympy.diff(dgdsxx, sxz)
ddgdsxxdsyz = sympy.diff(dgdsxx, syz)

ddgdsxxdsxx = ddgdsxxdsxx.subs([(F, 0.25216953733566727),
                                (G, 0.8254230293025175),
                                (H, 0.17457697069748246),
                                (L, 2.2380520016508463),
                                (M, 2.2380520016508463),
                                (N, 2.2380520016508463),
                                (sxx, 500), (syy, 700), (szz, 200),
                                (sxy, 100), (sxz, 200), (syz, 300), ])

ddgdsxxdsyy = ddgdsxxdsyy.subs([(F, 0.25216953733566727),
                                (G, 0.8254230293025175),
                                (H, 0.17457697069748246),
                                (L, 2.2380520016508463),
                                (M, 2.2380520016508463),
                                (N, 2.2380520016508463),
                                (sxx, 500), (syy, 700), (szz, 200),
                                (sxy, 100), (sxz, 200), (syz, 300), ])

ddgdsxxdszz = ddgdsxxdszz.subs([(F, 0.25216953733566727),
                                (G, 0.8254230293025175),
                                (H, 0.17457697069748246),
                                (L, 2.2380520016508463),
                                (M, 2.2380520016508463),
                                (N, 2.2380520016508463),
                                (sxx, 500), (syy, 700), (szz, 200),
                                (sxy, 100), (sxz, 200), (syz, 300), ])

ddgdsxxdsxy = ddgdsxxdsxy.subs([(F, 0.25216953733566727),
                                (G, 0.8254230293025175),
                                (H, 0.17457697069748246),
                                (L, 2.2380520016508463),
                                (M, 2.2380520016508463),
                                (N, 2.2380520016508463),
                                (sxx, 500), (syy, 700), (szz, 200),
                                (sxy, 100), (sxz, 200), (syz, 300), ])

ddgdsxxdsxz = ddgdsxxdsxz.subs([(F, 0.25216953733566727),
                                (G, 0.8254230293025175),
                                (H, 0.17457697069748246),
                                (L, 2.2380520016508463),
                                (M, 2.2380520016508463),
                                (N, 2.2380520016508463),
                                (sxx, 500), (syy, 700), (szz, 200),
                                (sxy, 100), (sxz, 200), (syz, 300), ])

ddgdsxxdsyz = ddgdsxxdsyz.subs([(F, 0.25216953733566727),
                                (G, 0.8254230293025175),
                                (H, 0.17457697069748246),
                                (L, 2.2380520016508463),
                                (M, 2.2380520016508463),
                                (N, 2.2380520016508463),
                                (sxx, 500), (syy, 700), (szz, 200),
                                (sxy, 100), (sxz, 200), (syz, 300), ])

print(ddgdsxxdsxx, ddgdsxxdsyy, ddgdsxxdszz, ddgdsxxdsxy, ddgdsxxdsxz, ddgdsxxdsyz)

# print(ddgdsxxdsxx.subs([(F, 0.25216953733566727),
#                         (G, 0.8254230293025175),
#                         (H, 0.17457697069748246),
#                         (L, 2.2380520016508463),
#                         (M, 2.2380520016508463),
#                         (N, 2.2380520016508463),
#                         (sxx, 500), (syy, 700), (szz, 200),
#                         (sxy, 100), (sxz, 200), (syz, 300), ]))
