#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <FPT.h>
#include "M3d_matrix_toolsS.c"

int window_size = 800;
int num_points[10000], num_polys[10000];
int psize[50][1000];
int con[100][1000][1000];
double x[10][50000], y[10][50000], z[10][50000];

int find_shape_center(int file_number) {
    // would just looping through all of the points in the polygon work as well/better? Or is this better? This can count points multiple times as they're used

    // find center through center of mass
    double x_total = 0, y_total = 0, z_total = 0;
    int count = 0;

    // get the x, y values of each point that is in each polygon
    for (int i = 0; i < num_polys[file_number]; i++) {
        for (int j = 0; j < psize[file_number][i]; j++) {
            x_total += x[file_number][con[file_number][i][j]];
            y_total += y[file_number][con[file_number][i][j]];
            z_total += z[file_number][con[file_number][i][j]];
            count   += 1;
        }
    }

    x[file_number][num_points[file_number]] = x_total / count;
    y[file_number][num_points[file_number]] = y_total / count;
    z[file_number][num_points[file_number]] = z_total / count;
}

void draw_polygon(int file_number) {
    // draw each polygon in the file
    for (int poly_num = 0; poly_num < num_polys[file_number]; poly_num++) {
        double x_points[1000], y_points[1000], z_points[1000];

        // get the x, y values of each point that is in each polygon
        for (int point = 0; point < psize[file_number][poly_num]; point++) {
            x_points[point] = x[file_number][con[file_number][poly_num][point]];
            y_points[point] = y[file_number][con[file_number][poly_num][point]];
            z_points[point] = z[file_number][con[file_number][poly_num][point]];
        }

        // do math to calculate 3d aspect of this shape
        double h;
        h = tan(60);

        for (int point = 0; point < psize[file_number][poly_num]; point++) {
            double x_bar, y_bar;
            x_bar = x_points[point] / z_points[point];
            y_bar = y_points[point] / z_points[point];

            double x_barbar, y_barbar;
            x_barbar = (window_size / 2) * (x_bar / h) + (window_size / 2);
            y_barbar = (window_size / 2) * (y_bar / h) + (window_size / 2);

            x_points[point] = x_barbar;
            y_points[point] = y_barbar;

        }

        // print polygon
        G_polygon(x_points, y_points, psize[file_number][poly_num]);
    }
}

int main(int argc, char **argv) {
	// get arguments on command line and check for files existence
    FILE *fp;
	int num_objects = argc - 1;

    // setup polygons, points, and numbers for each on the command line
   	for (int k = 0; k < num_objects; k++) {
		fp = fopen(argv[k + 1], "r");

        // if file doesn't exist, print error and quit
		if (fp == NULL) {
			printf("Can't open file, %s\n", argv[2]);
			exit(0);
		}

        // else, scan in values
        else {
		    // first line is the number of points
		    fscanf(fp, "%d", &num_points[k]);

            // next num_points lines are the x and y values of each point
		    for (int i = 0; i < num_points[k]; i++) {
                fscanf(fp, "%lf %lf %lf", &x[k][i], &y[k][i], &z[k][i]);
            }

            // next line is the number of polygons
		    fscanf(fp, "%d", &num_polys[k]);

            // next num_poly lines are the number of lines in this polygon
		    for (int i = 0; i < num_polys[k]; i++) {
			    fscanf(fp, "%d", &psize[k][i]);

                // followed by the points that are connected in this polygon
			    for (int j = 0; j < psize[k][i]; j++) {
                    fscanf(fp, "%d", &con[k][i][j]);
                }
		    }
        }
	}

    // make movement matrices
    // translate 1 in each direction matrices
    double x_trans_forw[4][4], y_trans_forw[4][4], z_trans_forw[4][4];
    double x_trans_back[4][4], y_trans_back[4][4], z_trans_back[4][4];
    M3d_make_translation(x_trans_forw,  0.1,  0,  0);
    M3d_make_translation(y_trans_forw,  0,  0.1,  0);
    M3d_make_translation(z_trans_forw,  0,  0,  0.1);
    M3d_make_translation(x_trans_back, -0.1,  0,  0);
    M3d_make_translation(y_trans_back,  0, -0.1,  0);
    M3d_make_translation(z_trans_back,  0,  0, -0.1);

    // rotate two degrees in each direction matrices
    double sine = sin((2 * M_PI) / 180);
    double cose = cos((2 * M_PI) / 180);
    double x_rot_forw[4][4], y_rot_forw[4][4], z_rot_forw[4][4];
    double x_rot_back[4][4], y_rot_back[4][4], z_rot_back[4][4];
    M3d_make_rotation_x(x_rot_forw,  sine,  cose);
    M3d_make_rotation_y(y_rot_forw,  sine,  cose);
    M3d_make_rotation_z(z_rot_forw,  sine,  cose);
    M3d_make_rotation_x(x_rot_back, -sine, cose);
    M3d_make_rotation_y(y_rot_back, -sine, cose);
    M3d_make_rotation_z(z_rot_back, -sine, cose);

	// setup graphics, draw black screen
    G_init_graphics(window_size, window_size);
    G_rgb(0, 0, 0);
    G_fill_rectangle(0, 0, window_size, window_size);

    // translate shapes to origin
    for (int file_number = 0; file_number < num_objects; file_number++) {
        // storing center in x, y, z, so that operations on the figure will also move the center
        find_shape_center(file_number);
        double translate_origin[4][4];

        M3d_make_translation(translate_origin, -x[file_number][num_points[file_number]], -y[file_number][num_points[file_number]], -z[file_number][num_points[file_number]]);
        M3d_mat_mult_points(x[file_number], y[file_number], z[file_number], translate_origin, x[file_number], y[file_number], z[file_number], num_points[file_number] + 1);
    }

    int file_number = 0;
    int last_num = 0;

    _Bool reverse   = False;
    _Bool rotate    = False;
    _Bool translate = False;

    draw_polygon(file_number);

    // loop for getting keys and drawing polygons
    while (1) {
        double transform[4][4];
        M3d_make_identity(transform);

        int key = G_wait_key();
        int key_num = key - 49;

        printf("Key_num is %d\n", key_num);

        // display new object
        if (key_num < num_objects) {
            file_number = key_num;
        }

        // reverse the transform/rotation
        else if (key == 'c') {
            if (reverse == True)    reverse = False;
            else                    reverse = True;
        }

        // press q to quit
        else if (key == 'q')        exit(0);

        // toggle rotating
        else if (key == 'r') {
            if (rotate == True) {
                rotate    = False;
                translate = True;
            }
            else {
                rotate    = True;
                translate = False;
            }
        }

        // toggle translating
        else if (key == 't') {
            if (translate == True) {
                rotate    = True;
                translate = False;
            }
            else {
                rotate    = False;
                translate = True;
            }
        }

        // make translation matrix and translate
        else if (translate == True) {
            if (key == 'x') {
                if      (reverse == True)   M3d_copy_mat(transform, x_trans_back);
                else if (reverse == False)  M3d_copy_mat(transform, x_trans_forw);
            }
            else if (key == 'y') {
                if      (reverse == True)   M3d_copy_mat(transform, y_trans_back);
                else if (reverse == False)  M3d_copy_mat(transform, y_trans_forw);
            }
            else if (key == 'z') {
                if      (reverse == True)   M3d_copy_mat(transform, z_trans_back);
                else if (reverse == False)  M3d_copy_mat(transform, z_trans_forw);
            }
        }

        // make rotation matrix
        else if (rotate == True) {
            // translate_center * rotate * translate_back
            M3d_make_translation(transform, -x[file_number][num_points[file_number]], -y[file_number][num_points[file_number]], -z[file_number][num_points[file_number]]);

            if (key == 'x') {
                if      (reverse == True)   M3d_mat_mult(transform, x_rot_back, transform);
                else if (reverse == False)  M3d_mat_mult(transform, x_rot_forw, transform);
            }
            else if (key == 'y') {
                if      (reverse == True)   M3d_mat_mult(transform, y_rot_back, transform);
                else if (reverse == False)  M3d_mat_mult(transform, y_rot_forw, transform);
            }
            else if (key == 'z') {
                if      (reverse == True)   M3d_mat_mult(transform, z_rot_back, transform);
                else if (reverse == False)  M3d_mat_mult(transform, z_rot_forw, transform);
            }

            // this should work as if it didnt fit other conditions it will continue anyways?
            else {
                continue;
            }

            // now translate back to original location?
            double translate_back[4][4];
            M3d_make_translation(translate_back, x[file_number][num_points[file_number]], y[file_number][num_points[file_number]], z[file_number][num_points[file_number]]);
            M3d_mat_mult(transform, translate_back, transform);
        }

        else {
            continue;
        }

        // multiply matrix by points
        M3d_mat_mult_points(x[file_number], y[file_number], z[file_number], transform, x[file_number], y[file_number], z[file_number], num_points[file_number] + 1);

        G_rgb(0, 0, 0);
	    G_fill_rectangle(0, 0, window_size, window_size);

        printf("Reverse %d\nRotate %d\nTranslate %d\n", reverse, rotate, translate);
        printf("Displaying %d\n\n", file_number);
        G_rgb(0, 1, 0);
        draw_polygon(file_number);
    }
}
