include("debug.jl")

xyz_to_gen_prop_idx_map = Dict(1 => 1, 2 => 5, 3 => 9)
gen_to_xyz_prop_idx_map = Dict(1 => 1, 5 => 2, 9 => 3)
dict_fd_to_str = Dict(
    1 => "fd.xx",
    2 => "fd.yy",
    3 => "fd.zz",
    4 => "fd.xy1",
    5 => "fd.xz1",
    6 => "fd.yz1",
    7 => "fd.yx1",
    8 => "fd.zx1",
    9 => "fd.zy1",
    10 => "fd.xy2",
    11 => "fd.xz2",
    12 => "fd.yz2",
    13 => "fd.yx2",
    14 => "fd.zx2",
    15 => "fd.zy2",
    16 => "fd.xy3",
    17 => "fd.xz3",
    18 => "fd.yz3",
    19 => "fd.yx3",
    20 => "fd.zx3",
    21 => "fd.zy3",
)


function loop()
    output = Dict()
    for n = 1:81
        xyz_n = XYZIndex_from_GenIndex(n)
        if xyz_n == -1
            continue
        end

        xyz_n_terms_P = Vector{String}()
        xyz_n_terms_PT = Vector{String}()

        d1, d2, d3, d4 = FlavorsFromGenIndex(n)
        for m = 1:81
            g1, g2, g3, g4 = FlavorsFromGenIndex(m)
            v13_flavor = GenIndexFromFlavors(d1, g2, d3, g3)
            v24_flavor = GenIndexFromFlavors(d2, g1, d4, g4)


            let v13_flavor_xyz = XYZIndex_from_GenIndex(v13_flavor),
                v24_flavor_xyz = XYZIndex_from_GenIndex(v24_flavor),
                P_idx_1 = g1 * 3 + g2 + 1,
                P_idx_2 = g3 * 3 + g4 + 1


                if (
                    v13_flavor_xyz != -1 &&
                    v24_flavor_xyz != -1 &&
                    P_idx_1 in values(xyz_to_gen_prop_idx_map) &&
                    P_idx_2 in values(xyz_to_gen_prop_idx_map)
                )
                    P_xyz_idx_1 = gen_to_xyz_prop_idx_map[P_idx_1]
                    P_xyz_idx_2 = gen_to_xyz_prop_idx_map[P_idx_2]
                    push!(
                        xyz_n_terms_P,
                        "V13[$(dict_fd_to_str[v13_flavor_xyz])] * V24[$(dict_fd_to_str[v24_flavor_xyz])] * P_($P_xyz_idx_1,$P_xyz_idx_2)",
                    )
                end
            end
        end

        for m = 1:81
            g1, g2, g3, g4 = FlavorsFromGenIndex(m)


            v31_flavor = GenIndexFromFlavors(d1, g3, d3, g2)
            v42_flavor = GenIndexFromFlavors(d2, g4, d4, g1)

            PT_idx_1 = g1 * 3 + g2 + 1
            PT_idx_2 = g3 * 3 + g4 + 1

            let v31_flavor_xyz = XYZIndex_from_GenIndex(v31_flavor),
                v42_flavor_xyz = XYZIndex_from_GenIndex(v42_flavor),
                PT_idx_1 = g1 * 3 + g2 + 1,
                PT_idx_2 = g3 * 3 + g4 + 1

                if (
                    v31_flavor_xyz != -1 &&
                    v42_flavor_xyz != -1 &&
                    PT_idx_1 in values(xyz_to_gen_prop_idx_map) &&
                    PT_idx_2 in values(xyz_to_gen_prop_idx_map)
                )

                    PT_xyz_idx_1 = gen_to_xyz_prop_idx_map[PT_idx_1]
                    PT_xyz_idx_2 = gen_to_xyz_prop_idx_map[PT_idx_2]

                    push!(
                        xyz_n_terms_PT,
                        "V31[$(dict_fd_to_str[v31_flavor_xyz])] * V42[$(dict_fd_to_str[v42_flavor_xyz])] * PT_($PT_xyz_idx_1,$PT_xyz_idx_2)",
                    )
                end
            end
        end
        output[xyz_n] = (xyz_n_terms_P, xyz_n_terms_PT)
    end
    return output
end

function print_loop(loop_output)
    for xyz_n = 1:21
        flavour_string = dict_fd_to_str[xyz_n]
        print("X_sum[21+$flavour_string] += (")
        (xyz_n_terms_P, xyz_n_terms_PT) = loop_output[xyz_n]
        print("(")
        println(join(xyz_n_terms_P, " +\n\t"))
        print(")")
        println("+")
        print("(")
        println(join(xyz_n_terms_PT, " +\n\t"))
        print(")")
        println(")")
    end

end


print_loop(loop())
