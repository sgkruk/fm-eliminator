#Fourier-Motzkin elimination

# Introduction #

A user-friendly FM eliminator.


# Details #

The idea is to take easy to read linear inequalities, a set of variables on which to project and produce an easy to read projection.

For example if the input file is
```
y11+y12+y13 = 1
y21+y22+y23 = 1

y11+y21 <=1
y12+y22 <=1
y13+y23 <=1

y11>=0
y12>=0
y13>=0
y21>=0
y22>=0
y23>=0

1y11+2y12+3y13 - x1 = 0
1y21+2y22+3y23 - x2 = 0
```
and that you ask to project onto the x-space (x1,x2), you obtain

```
-1.00000 x1   -1.00000 x2   <=    -3
              -1.00000 x2   <=    -1
+1.00000 x1   +1.00000 x2   <=     5
+1.00000 x1                 <=     3
-1.00000 x1                 <=    -1
              +1.00000 x2   <=     3
```