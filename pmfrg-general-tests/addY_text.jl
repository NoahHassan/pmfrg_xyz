X_sum[21+fd.xx] += (
    (
        V13[fd.xx] * V24[fd.xx] * P_(1, 1) +
        V13[fd.xy2] * V24[fd.xy2] * P_(2, 2) +
        V13[fd.xz2] * V24[fd.xz2] * P_(3, 3)
    ) + (
        V31[fd.xx] * V42[fd.xx] * PT_(1, 1) +
        V31[fd.xy2] * V42[fd.xy2] * PT_(2, 2) +
        V31[fd.xz2] * V42[fd.xz2] * PT_(3, 3)
    )
)

X_sum[21+fd.yy] += (
    (
        V13[fd.yy] * V24[fd.yy] * P_(2, 2) +
        V13[fd.yx2] * V24[fd.yx2] * P_(1, 1) +
        V13[fd.yz2] * V24[fd.yz2] * P_(3, 3)
    ) + (
        V31[fd.yy] * V42[fd.yy] * PT_(2, 2) +
        V31[fd.yx2] * V42[fd.yx2] * PT_(1, 1) +
        V31[fd.yz2] * V42[fd.yz2] * PT_(3, 3)
    )
)

X_sum[21+fd.zz] += (
    (
        V13[fd.zz] * V24[fd.zz] * P_(3, 3) +
        V13[fd.zx2] * V24[fd.zx2] * P_(1, 1) +
        V13[fd.zy2] * V24[fd.zy2] * P_(2, 2)
    ) + (
        V31[fd.zz] * V42[fd.zz] * PT_(3, 3) +
        V31[fd.zx2] * V42[fd.zx2] * PT_(1, 1) +
        V31[fd.zy2] * V42[fd.zy2] * PT_(2, 2)
    )
)

### Yab1 = Vab3 Vab3 + Vab1 Vab1 + (w -- -w + t)

X_sum[21+fd.xy1] += (
    (V13[fd.xy3] * V24[fd.xy3] * P_(2, 1) + V13[fd.xy1] * V24[fd.xy1] * P_(1, 2)) +
    (V31[fd.xy3] * V42[fd.xy3] * PT_(1, 2) + V31[fd.xy1] * V42[fd.xy1] * PT_(2, 1))
)

X_sum[21+fd.xz1] += (
    (V13[fd.xz3] * V24[fd.xz3] * P_(3, 1) + V13[fd.xz1] * V24[fd.xz1] * P_(1, 3)) +
    (V31[fd.xz3] * V42[fd.xz3] * PT_(1, 3) + V31[fd.xz1] * V42[fd.xz1] * PT_(3, 1))
)

X_sum[21+fd.yx1] += (
    (V13[fd.yx3] * V24[fd.yx3] * P_(1, 2) + V13[fd.yx1] * V24[fd.yx1] * P_(2, 1)) +
    (V31[fd.yx3] * V42[fd.yx3] * PT_(2, 1) + V31[fd.yx1] * V42[fd.yx1] * PT_(1, 2))
)

X_sum[21+fd.yz1] += (
    (V13[fd.yz3] * V24[fd.yz3] * P_(3, 2) + V13[fd.yz1] * V24[fd.yz1] * P_(2, 3)) +
    (V31[fd.yz3] * V42[fd.yz3] * PT_(2, 3) + V31[fd.yz1] * V42[fd.yz1] * PT_(3, 2))
)

X_sum[21+fd.zx1] += (
    (V13[fd.zx3] * V24[fd.zx3] * P_(1, 3) + V13[fd.zx1] * V24[fd.zx1] * P_(3, 1)) +
    (V31[fd.zx3] * V42[fd.zx3] * PT_(3, 1) + V31[fd.zx1] * V42[fd.zx1] * PT_(1, 3))
)

X_sum[21+fd.zy1] += (
    (V13[fd.zy3] * V24[fd.zy3] * P_(2, 3) + V13[fd.zy1] * V24[fd.zy1] * P_(3, 2)) +
    (V31[fd.zy3] * V42[fd.zy3] * PT_(3, 2) + V31[fd.zy1] * V42[fd.zy1] * PT_(2, 3))
)

### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

X_sum[21+fd.xy2] += (
    (
        V13[fd.xx] * V24[fd.yx2] * P_(1, 1) +
        V13[fd.xy2] * V24[fd.yy] * P_(2, 2) +
        V13[fd.xz2] * V24[fd.yz2] * P_(3, 3)
    ) + (
        V31[fd.xx] * V42[fd.yx2] * PT_(1, 1) +
        V31[fd.xy2] * V42[fd.yy] * PT_(2, 2) +
        V31[fd.xz2] * V42[fd.yz2] * PT_(3, 3)
    )
)

X_sum[21+fd.xz2] += (
    (
        V13[fd.xx] * V24[fd.zx2] * P_(1, 1) +
        V13[fd.xz2] * V24[fd.zz] * P_(3, 3) +
        V13[fd.xy2] * V24[fd.zy2] * P_(2, 2)
    ) + (
        V31[fd.xx] * V42[fd.zx2] * PT_(1, 1) +
        V31[fd.xz2] * V42[fd.zz] * PT_(3, 3) +
        V31[fd.xy2] * V42[fd.zy2] * PT_(2, 2)
    )
)

X_sum[21+fd.yx2] += (
    (
        V13[fd.yy] * V24[fd.xy2] * P_(2, 2) +
        V13[fd.yx2] * V24[fd.xx] * P_(1, 1) +
        V13[fd.yz2] * V24[fd.xz2] * P_(3, 3)
    ) + (
        V31[fd.yy] * V42[fd.xy2] * PT_(2, 2) +
        V31[fd.yx2] * V42[fd.xx] * PT_(1, 1) +
        V31[fd.yz2] * V42[fd.xz2] * PT_(3, 3)
    )
)

X_sum[21+fd.yz2] += (
    (
        V13[fd.yy] * V24[fd.zy2] * P_(2, 2) +
        V13[fd.yz2] * V24[fd.zz] * P_(3, 3) +
        V13[fd.yx2] * V24[fd.zx2] * P_(1, 1)
    ) + (
        V31[fd.yy] * V42[fd.zy2] * PT_(2, 2) +
        V31[fd.yz2] * V42[fd.zz] * PT_(3, 3) +
        V31[fd.yx2] * V42[fd.zx2] * PT_(1, 1)
    )
)

X_sum[21+fd.zx2] += (
    (
        V13[fd.zz] * V24[fd.xz2] * P_(3, 3) +
        V13[fd.zx2] * V24[fd.xx] * P_(1, 1) +
        V13[fd.zy2] * V24[fd.xy2] * P_(2, 2)
    ) + (
        V31[fd.zz] * V42[fd.xz2] * PT_(3, 3) +
        V31[fd.zx2] * V42[fd.xx] * PT_(1, 1) +
        V31[fd.zy2] * V42[fd.xy2] * PT_(2, 2)
    )
)

X_sum[21+fd.zy2] += (
    (
        V13[fd.zz] * V24[fd.yz2] * P_(3, 3) +
        V13[fd.zy2] * V24[fd.yy] * P_(2, 2) +
        V13[fd.zx2] * V24[fd.yx2] * P_(1, 1)
    ) + (
        V31[fd.zz] * V42[fd.yz2] * PT_(3, 3) +
        V31[fd.zy2] * V42[fd.yy] * PT_(2, 2) +
        V31[fd.zx2] * V42[fd.yx2] * PT_(1, 1)
    )
)

### Yab3 = Vab3 Vba1 + Vab1 Vba3 + (w -- -w + t)

X_sum[21+fd.xy3] += (
    (V13[fd.xy3] * V24[fd.yx1] * P_(2, 1) + V13[fd.xy1] * V24[fd.yx3] * P_(1, 2)) +
    (V31[fd.xy3] * V42[fd.yx1] * PT_(1, 2) + V31[fd.xy1] * V42[fd.yx3] * PT_(2, 1))
)

X_sum[21+fd.xz3] += (
    (V13[fd.xz3] * V24[fd.zx1] * P_(3, 1) + V13[fd.xz1] * V24[fd.zx3] * P_(1, 3)) +
    (V31[fd.xz3] * V42[fd.zx1] * PT_(1, 3) + V31[fd.xz1] * V42[fd.zx3] * PT_(3, 1))
)
X_sum[21+fd.yx3] += (
    (V13[fd.yx3] * V24[fd.xy1] * P_(1, 2) + V13[fd.yx1] * V24[fd.xy3] * P_(2, 1)) +
    (V31[fd.yx3] * V42[fd.xy1] * PT_(2, 1) + V31[fd.yx1] * V42[fd.xy3] * PT_(1, 2))
)
X_sum[21+fd.yz3] += (
    (V13[fd.yz3] * V24[fd.zy1] * P_(3, 2) + V13[fd.yz1] * V24[fd.zy3] * P_(2, 3)) +
    (V31[fd.yz3] * V42[fd.zy1] * PT_(2, 3) + V31[fd.yz1] * V42[fd.zy3] * PT_(3, 2))
)

X_sum[21+fd.zx3] += (
    (V13[fd.zx3] * V24[fd.xz1] * P_(1, 3) + V13[fd.zx1] * V24[fd.xz3] * P_(3, 1)) +
    (V31[fd.zx3] * V42[fd.xz1] * PT_(3, 1) + V31[fd.zx1] * V42[fd.xz3] * PT_(1, 3))
)

X_sum[21+fd.zy3] += (
    (V13[fd.zy3] * V24[fd.yz1] * P_(2, 3) + V13[fd.zy1] * V24[fd.yz3] * P_(3, 2)) +
    (V31[fd.zy3] * V42[fd.yz1] * PT_(3, 2) + V31[fd.zy1] * V42[fd.yz3] * PT_(2, 3))
)
