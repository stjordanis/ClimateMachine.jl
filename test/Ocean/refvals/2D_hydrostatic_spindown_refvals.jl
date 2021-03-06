# [
#  [ MPIStateArray Name, Field Name, Maximum, Minimum, Mean, Standard Deviation ],
#  [         :                :          :        :      :          :           ],
# ]

parr = [
    ["2D state", "η", 12, 12, 0, 12],
    ["2D state", "U[1]", 12, 12, 0, 12],
    ["2D state", "U[2]", 0, 0, 0, 0],
]

explicit = [
    [
        "2D state",
        "η",
        -8.52722969951589915283e-01,
        8.52846676313531282254e-01,
        -2.49578135935735214742e-16,
        6.03454239990563690021e-01,
    ],
    [
        "2D state",
        "U[1]",
        -3.15431401945821825450e+01,
        3.15431401945818628008e+01,
        6.11504145930918957291e-15,
        2.24273815174625497093e+01,
    ],
    [
        "2D state",
        "U[2]",
        -7.62224398365580242501e-13,
        9.72156930292624284356e-13,
        1.39269607441935025982e-14,
        1.95606703846656748360e-13,
    ],
]

refVals = (explicit = (explicit, parr),)
