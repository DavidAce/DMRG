# With $\chi_\text{max} = 16$ for 4000 iDMRG steps  

## ncv = 8

| Method                             |  iter  | time per step | total time |
| ---------------------------------- | :----: | ------------: | ---------: |
| Original (L theta HA HB R)         |   63   |          0.07 |      166.8 |
| mytens (HA HB L R)  idx({1},{0})   |   76   |          0.38 |       1400 |
| mytens (HA HB L R)  idx({0},{0})   | ~50-60 |             1 |      ~4000 |
| mymat (rowmajor) mytens(HA HB L R) |   68   |          0.25 |       1000 |
| mymat (colmajor) mytens(HA HB L R) |   70   |          0.32 |       1280 |
| mymat(rowmajor) (L HA HB R)        |   68   |          0.26 |       1210 |
| Original (theta L R HA HB)         |   70   |          0.19 |        700 |
| Original(theta R HB HA L)          |   58   |          0.07 |        163 |
| Original(theta L HA HB R)          |   58   |          0.75 |        169 |





Number of contractions for each ncv



| ncv  | lanczos iter | contraction counter | ratio | time per step |
| ---- | ------------ | ------------------- | ----- | ------------- |
| 4    | 594          | 1191                | ~2    | 0.33          |
| 6    | 121          | 367                 | ~3    | 0.01          |
| 8    | 62           | 253                 | ~4    | 0.07          |
| 10   | 36           | 186                 | ~5    | 0.057         |
| 12   | 26           | 163                 | ~6.27 | 0.52          |
| 14   |              |                     |       |               |
|      |              |                     |       |               |
|      |              |                     |       |               |
|      |              |                     |       |               |

