#include <stdio.h>
#include <math.h>

// here's how to print a 3x3 array
int M3d_print_mat (double a[4][4]) {
  for (int r = 0; r < 4; r++) {
      for (int c = 0; c < 4; c++) {
           printf(" %12.4lf ", a[r][c]);
      }
      printf("\n");
  }

  return 1;
}

// copy matrix
int M3d_copy_mat (double a[4][4], double b[4][4]) {
    // a = b
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            a[i][j] = b[i][j];
        }
    }

    return 1;
}

// all spots are 0, unless row = column
int M3d_make_identity (double a[4][4]) {
    // a = I
    int r, c;
    for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++) {
            if (r == c) {
                a[r][c] = 1.0;
            }
            else {
                a[r][c] = 0.0;
            }
        }
    }

    return 1;
}

// do a translation
int M3d_make_translation (double a[4][4], double dx, double dy, double dz) {
    M3d_make_identity(a) ;

    a[0][3] = dx;
    a[1][3] = dy;
    a[2][3] = dz;

    return 1;
}

// multiply two 3x3 matrices
int M3d_mat_mult (double res[4][4], double a[4][4], double b[4][4]) {
    // res = a * b
    // this is SAFE, i.e. the user can make a call such as
    // M3d_mat_mult(p, p, q) or M3d_mat_mult(p, q, p) or  M3d_mat_mult(p, p, p)
    double first[4][4];
    double second[4][4];

    M3d_copy_mat(first, a);
    M3d_copy_mat(second,b);


    // math
    for (int first_row = 0; first_row < 4; first_row++) {
        for (int second_col = 0; second_col < 4; second_col++) {
            double result_point = 0;

            for (int point = 0; point < 4; point++) {
                result_point += first[first_row][point] * second[point][second_col];
            }

            res[first_row][second_col] = result_point;
        }
    }

    return 1;
}

// make a scale
int M3d_make_scaling (double a[4][4], double sx, double sy, double sz) {
    M3d_make_identity(a);
    a[0][0] = sx;
    a[1][1] = sy;
    a[2][2] = sz;

    return 1;
}

// make rotation in radians
int M3d_make_rotation_x (double a[4][4], double sine, double cosine) {
    M3d_make_identity(a);

    // [1      0      0      0]
    // [0    cos   -sin      0]
    // [0    sin    cos      0]
    // [0      0      0      1]

    a[0][0] = 1;
    a[0][1] = 0;
    a[0][2] = 0;
    a[0][3] = 0;
    a[1][0] = 0;
    a[1][1] = cosine;
    a[1][2] = -sine;
    a[1][3] = 0;
    a[2][0] = 0;
    a[2][1] = sine;
    a[2][2] = cosine;
    a[2][3] = 0;
    a[3][0] = 0;
    a[3][1] = 0;
    a[3][2] = 0;
    a[3][3] = 1;

    return 1;
}

int M3d_make_rotation_y (double a[4][4], double sine, double cosine) {
    M3d_make_identity(a);

    // [ cos    0   sin     0]
    // [   0    1     0     0]
    // [-sin    0   cos     0]
    // [   0    0     0     1]

    a[0][0] = cosine;
    a[0][1] = 0;
    a[0][2] = sine;
    a[0][3] = 0;
    a[1][0] = 0;
    a[1][1] = 1;
    a[1][2] = 0;
    a[1][3] = 0;
    a[2][0] = -sine;
    a[2][1] = 0;
    a[2][2] = cosine;
    a[2][3] = 0;
    a[3][0] = 0;
    a[3][1] = 0;
    a[3][2] = 0;
    a[3][3] = 1;

    return 1;
}

int M3d_make_rotation_z (double a[4][4], double sine, double cosine) {
    M3d_make_identity(a);

    // [cos  -sin    0    0]
    // [sin   cos    0    0]
    // [  0     0    1    0]
    // [  0     0    0    1]

    a[0][0] = cosine;
    a[0][1] = -sine;
    a[0][2] = 0;
    a[0][3] = 0;
    a[1][0] = sine;
    a[1][1] = cosine;
    a[1][2] = 0;
    a[1][3] = 0;
    a[2][0] = 0;
    a[2][1] = 0;
    a[2][2] = 1;
    a[2][3] = 0;
    a[3][0] = 0;
    a[3][1] = 0;
    a[3][2] = 0;
    a[3][3] = 1;

    return 1;
}

int M3d_mat_mult_pt (double P[3], double m[4][4], double Q[3]) {
    // P = m*Q
    // SAFE, user may make a call like M3d_mat_mult_pt (W, m,W) ;

    // assume that the extra spot is a 1
    // [ 1,  2,  3]   [13]
    // [ 4,  5,  6] * [14]
    // [ 7,  8,  9]   [15]
    // [10, 11, 12]   [16]

    double point[3];
    point[0] = Q[0];
    point[1] = Q[1];
    point[2] = Q[2];

    P[0] = m[0][0] * point[0] + m[0][1] * point[1] + m[0][2] * point[2] + m[0][3];
    P[1] = m[1][0] * point[0] + m[1][1] * point[1] + m[1][2] * point[2] + m[1][3];
    P[2] = m[2][0] * point[0] + m[2][1] * point[1] + m[2][2] * point[2] + m[2][3];

    return 1;
}

int M3d_mat_mult_points (double *X, double *Y, double *Z, double m[4][4], double *x, double *y, double *z, int numpoints) {
    // |X0 X1 X2 ...|       |x0 x1 x2 ...|
    // |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
    // |Z0 Z1 Z2 ...|       |z0 z1 z2 ...|
    // | 1  1  1 ...|       | 1  1  1 ...|
    // SAFE, user may make a call like M3d_mat_mult_points (x,y,z, m, x,y,z, n) ;
    double t, u, v;

    for (int i = 0; i < numpoints; i++) {
        t = m[0][0]*x[i] + m[0][1]*y[i] + m[0][2]*z[i] + m[0][3];
        u = m[1][0]*x[i] + m[1][1]*y[i] + m[1][2]*z[i] + m[1][3];
        v = m[2][0]*x[i] + m[2][1]*y[i] + m[2][2]*z[i] + m[2][3];

        X[i] = t;
        Y[i] = u;
        Z[i] = v;
    }

    return 1;
}

/*int main() {
    double a[3][3], b[3][3], result[3][3];
    a[0][0] = 1;
    a[0][1] = 2;
    a[0][2] = 3;
    a[1][0] = 4;
    a[1][1] = 5;
    a[1][2] = 6;
    a[2][0] = 1;
    a[2][1] = 1;
    a[2][2] = 1;

    b[0][0] = 7;
    b[0][1] = 8;
    b[0][2] = 2;
    b[1][0] = 9;
    b[1][1] = 10;
    b[1][2] = 2;
    b[2][0] = 11;
    b[2][1] = 12;
    b[2][2] = 2;

    M3d_mat_mult(result, a, b);

    M3d_print_mat(result);

    printf("\n");

    M3d_make_rotation_cs(a, 1, 0);

    M3d_print_mat(a);
}*/
