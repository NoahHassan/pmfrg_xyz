using JLD2

### Step by step

### Step one: First iteration, after getXBubble:

noahX = load_object("noahX.jld2")
yannX = load_object("yannikX.jld2")

### Checks X_a = X_a
for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahX[n, Rij, is, it, iu], yannX.x[1][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahX[1, Rij, is, it, iu])
                println(yannX.x[1][Rij, is, it, iu])
            end
        end
    end
end

### Checks X_b = X_b
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahX[3 + n, Rij, is, it, iu], yannX.x[2][Rij, is, it, iu], atol=1e-20))
                println("Conflict at n=$n, ($Rij, $is, $it, $iu)")
                println(noahX[3 + n, Rij, is, it, iu])
                println(yannX.x[2][Rij, is, it, iu])
            end
        end
    end
end

### Checks X_c = X_c
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahX[9 + n, Rij, is, it, iu], yannX.x[3][Rij, is, it, iu], atol=1e-20))
                println("Conflict at n=$n, ($Rij, $is, $it, $iu)")
                println(noahX[9 + n, Rij, is, it, iu])
                println(yannX.x[3][Rij, is, it, iu])
            end
        end
    end
end

### Checks Y_aa(s, t, u) = -X_ta(t, s, u)
for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahX[21 + n, Rij, is, it, iu], -yannX.x[4][Rij, it, is, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahX[21 + n, Rij, is, it, iu])
                println(yannX.x[4][Rij, it, is, iu])
            end
        end
    end
end

### Checks Y_ab1(s, t, u) = -X_tc(t, s, u)
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahX[21 + 3 + n, Rij, is, it, iu], -yannX.x[6][Rij, it, is, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahX[21 + 3 + n, Rij, is, it, iu])
                println(yannX.x[6][Rij, it, is, iu])
            end
        end
    end
end

### Checks Y_ab2(s, t, u) = -X_tb(t, s, u)
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahX[21 + 9 + n, Rij, is, it, iu], -yannX.x[5][Rij, it, is, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahX[21 + 9 + n, Rij, is, it, iu])
                println(yannX.x[5][Rij, it, is, iu])
            end
        end
    end
end

############# SUCCESS #############

### Step two: First iteration, after symmetrizeBubble

noahXSymm = load_object("noahXSymm.jld2")
yannXSymm = load_object("yannikXSymm.jld2")

### Checks X_a = X_a
for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahXSymm[n, Rij, is, it, iu], yannXSymm.x[1][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahXSymm[1, Rij, is, it, iu])
                println(yannXSymm.x[1][Rij, is, it, iu])
            end
        end
    end
end

### Checks X_b = X_b
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahXSymm[3 + n, Rij, is, it, iu], yannXSymm.x[2][Rij, is, it, iu], atol=1e-20))
                println("Conflict at n=$n, ($Rij, $is, $it, $iu)")
                println(noahXSymm[3 + n, Rij, is, it, iu])
                println(yannXSymm.x[2][Rij, is, it, iu])
            end
        end
    end
end

### Checks X_c = X_c
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahXSymm[9 + n, Rij, is, it, iu], yannXSymm.x[3][Rij, is, it, iu], atol=1e-20))
                println("Conflict at n=$n, ($Rij, $is, $it, $iu)")
                println(noahXSymm[9 + n, Rij, is, it, iu])
                println(yannXSymm.x[3][Rij, is, it, iu])
            end
        end
    end
end

### Checks Y_aa(s, t, u) = -X_ta(t, s, u)
for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahXSymm[21 + n, Rij, is, it, iu], -yannXSymm.x[4][Rij, it, is, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahXSymm[21 + n, Rij, is, it, iu])
                println(yannXSymm.x[4][Rij, it, is, iu])
            end
        end
    end
end

### Checks Y_ab1(s, t, u) = -X_tc(t, s, u)
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahXSymm[21 + 3 + n, Rij, is, it, iu], -yannXSymm.x[6][Rij, it, is, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahXSymm[21 + 3 + n, Rij, is, it, iu])
                println(yannXSymm.x[6][Rij, it, is, iu])
            end
        end
    end
end

### Checks Y_ab2(s, t, u) = -X_tb(t, s, u)
for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahXSymm[21 + 9 + n, Rij, is, it, iu], -yannXSymm.x[5][Rij, it, is, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahXSymm[21 + 9 + n, Rij, is, it, iu])
                println(yannXSymm.x[5][Rij, it, is, iu])
            end
        end
    end
end

############# SUCCESS #############

### Step three: First iteration, after addToVertexFromBubble

noahDer = load_object("noahDeriv.jld2")
yannDer = load_object("yannikDeriv.jld2")

for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahDer[n, Rij, is, it, iu], yannDer.x[1][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahDer[n, Rij, is, it, iu])
                println(yannDer.x[1][Rij, is, it, iu])
            end
        end
    end
end

for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahDer[3 + n, Rij, is, it, iu], yannDer.x[2][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahDer[3 + n, Rij, is, it, iu])
                println(yannDer.x[2][Rij, is, it, iu])
            end
        end
    end
end

for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahDer[9 + n, Rij, is, it, iu], yannDer.x[3][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahDer[9 + n, Rij, is, it, iu])
                println(yannDer.x[3][Rij, is, it, iu])
            end
        end
    end
end

############# SUCCESS #############

### Step four: First iteration, after symmetrizeVertex

noahDerSymm = load_object("noahDerivSymm.jld2")
yannDerSymm = load_object("yannikDerivSymm.jld2")

for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahDerSymm[n, Rij, is, it, iu], yannDerSymm.x[1][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahDerSymm[n, Rij, is, it, iu])
                println(yannDerSymm.x[1][Rij, is, it, iu])
            end
        end
    end
end

for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahDerSymm[3 + n, Rij, is, it, iu], yannDerSymm.x[2][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahDerSymm[3 + n, Rij, is, it, iu])
                println(yannDerSymm.x[2][Rij, is, it, iu])
            end
        end
    end
end

for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahDerSymm[9 + n, Rij, is, it, iu], yannDerSymm.x[3][Rij, is, it, iu], atol=1e-20))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahDerSymm[9 + n, Rij, is, it, iu])
                println(yannDerSymm.x[3][Rij, is, it, iu])
            end
        end
    end
end

############# SUCCESS #############

### Step five: Second iteration, initial Gamma

noahGam = load_object("noahGammaIt2.jld2")
yannGam = load_object("yannikGammaIt2.jld2")

for n in 1:3
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahGam[n, Rij, is, it, iu], yannGam.x[1][Rij, is, it, iu], atol=1e-15))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahGam[n, Rij, is, it, iu])
                println(yannGam.x[1][Rij, is, it, iu])
            end
        end
    end
end

for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahGam[3 + n, Rij, is, it, iu], yannGam.x[2][Rij, is, it, iu], atol=1e-15))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahGam[3 + n, Rij, is, it, iu])
                println(yannGam.x[2][Rij, is, it, iu])
            end
        end
    end
end

for n in 1:6
    for Rij in 1:2
        for is in 1:8, it in 1:8, iu in 1:8
            if(!isapprox(noahGam[9 + n, Rij, is, it, iu], yannGam.x[3][Rij, is, it, iu], atol=1e-15))
                println("Conflict at ($Rij, $is, $it, $iu)")
                println(noahGam[9 + n, Rij, is, it, iu])
                println(yannGam.x[3][Rij, is, it, iu])
            end
        end
    end
end