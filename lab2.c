#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <FPT.h>
#include "M2d_matrix_toolsS.c"

/**
Criticism:
    Currently, it calculates the center every time, and the matrices every time
    Instead, calculate the center once, and store that somewhere
    For every rotation, simply multiply the rotation matrix by the matrix the matrix that is the result of the two translation matrices and scale matrices
    Have that matrix stored somewhere, so that you don't need to build these matrices everytime
*/

int num_points[100];
int num_polys[100];
int psize[50][100];
int con[100][100][100];
double scale[10];
double center[10][2];
double x[10][100];
double y[10][100];
double red[10][100];
double grn[10][100];
double blu[10][100];

int find_scale(int file_number) {
    // find scale by finding the bounds
    double min_x = x[file_number][con[file_number][0][0]];
    double min_y = y[file_number][con[file_number][0][0]];
    double max_x = x[file_number][con[file_number][0][0]];
    double max_y = y[file_number][con[file_number][0][0]];

    for (int i = 0; i < num_polys[file_number]; i++) {
        double poly_x[10];
        double poly_y[10];

        // get the x, y values of each point that is in each polygon
        for (int j = 0; j < psize[file_number][i]; j++) {
            poly_x[j] = x[file_number][con[file_number][i][j]];
            poly_y[j] = y[file_number][con[file_number][i][j]];

            if (poly_x[j] < min_x) {
                min_x = poly_x[j];
            }
            else if (poly_x[j] > max_x) {
                max_x = poly_x[j];
            }
            if (poly_y[j] < min_y) {
                min_y = poly_y[j];
            }
            else if (poly_y[j] > max_y) {
                max_y = poly_y[j];
            }
        }
    }

    scale[file_number] = 300 / (max_x - center[file_number][0]);
    if ((300 / (max_y - center[file_number][1])) < scale[file_number]) {
        scale[file_number] = 300 / (max_y - center[file_number][1]);
    }

    printf("found scale\n");

    return 1;
}

int find_center(int file_number) {
    // find center through average
    // center of mass
    double x_avg = 0;
    double y_avg = 0;
    int count = 0;

    // get the x, y values of each point that is in each polygon
    for (int i = 0; i < num_polys[file_number]; i++) {
        for (int j = 0; j < psize[file_number][i]; j++) {
            x_avg += x[file_number][con[file_number][i][j]];
            y_avg += y[file_number][con[file_number][i][j]];
            count += 1;
        }
    }

    center[file_number][0] = x_avg / count;
    center[file_number][1] = y_avg / count;

    printf("center is %lf, %lf\n", center[file_number][0], center[file_number][1]);

    return 1;
}

void all_moves(double poly_points[2][100], int length, int file_number, double cosine, double sine) {
    // makes matrix to translate polygon to origin
    double translate_origin[3][3];
    M2d_make_translation(translate_origin, -center[file_number][0], -center[file_number][1]);

    // makes matrix to scale polygon to size of screen
    double scaled[3][3];
    M2d_make_scaling(scaled, scale[file_number], scale[file_number]);

    // makes matrix to rotate polygon to desired rotation
    double rotation[3][3];
    M2d_make_rotation_cs(rotation, cosine, sine);

    // makes matrix to translate polygon to center
    double translation_center[3][3];
    M2d_make_translation(translation_center, 300, 300);

    // multiplies all matrices together to get final result
    double result[3][3];
    M2d_mat_mult(result, scaled, translate_origin);
    M2d_mat_mult(result, rotation, result);
    M2d_mat_mult(result, translation_center, result);
    M2d_mat_mult_points(poly_points[0], poly_points[1], result, poly_points[0], poly_points[1], length);
}

void draw_polygon(int file_number, double radians) {
    double sine = sin(radians);
    double cosine = cos(radians);

    // draw each polygon in each file
    for (int poly_num = 0; poly_num < num_polys[file_number]; poly_num++) {
        double poly_points[2][100];
        int length = psize[file_number][poly_num];

        // get the x, y values of each point that is in each polygon
        for (int point = 0; point < length; point++) {
            poly_points[0][point] = x[file_number][con[file_number][poly_num][point]];
            poly_points[1][point] = y[file_number][con[file_number][poly_num][point]];
        }

        all_moves(poly_points, length, file_number, cosine, sine);

        // set RGB and print polygon
        G_rgb(red[file_number][poly_num], grn[file_number][poly_num], blu[file_number][poly_num]);
        G_fill_polygon(poly_points[0], poly_points[1], psize[file_number][poly_num]);
    }
}

void setup_polygons(int num_objects) {
    for (int file_number = 0; file_number < num_objects; file_number++) {
        find_center(file_number);
        find_scale(file_number);
    }
    printf("centered and scaled\n");
}

int main(int argc, char **argv) {
	// get arguments on command line and check for files existence
    FILE *fp;
	int num_objects = argc - 1;

   	for (int k = 0; k < num_objects; k++) {
		fp = fopen(argv[k + 1], "r");

		if (fp == NULL) {
			printf("Can't open file, %s\n", argv[2]);
			exit(0);
		}
        // first line is the number of points
        fscanf(fp, "%d", &num_points[k]);

        //printf("processed num_points\n");

        // next NUMPOINTS lines are the x and y values of each point
		for (int i = 0; i < num_points[k]; i++) {
			fscanf(fp, "%lf %lf", &x[k][i], &y[k][i]);
		}

        //printf("processed points\n");

        // next line is the number of polygons
		fscanf(fp, "%d", &num_polys[k]);

        //printf("processed num_polys\n");

        // next NUMPOLY lines are the number of lines in this polygon
		for (int i = 0; i < num_polys[k]; i++) {
			fscanf(fp, "%d", &psize[k][i]);

            // followed by the points that are connected in this polygon
			for (int j = 0; j < psize[k][i]; j++) {
				fscanf(fp, "%d", &con[k][i][j]);
			}
		}

        //printf("processed lines\n");

        // next NUMPOLY lines are the RGB values of each polygon
		for (int i = 0; i < num_polys[k]; i++) {
			fscanf(fp, "%lf %lf %lf", &red[k][i], &grn[k][i], &blu[k][i]);
		}

        //printf("processed RGB\n");
	}

	// draw these polygons!
    G_init_graphics(600, 600);
    double radians[num_objects];

    for (int num_object = 0; num_object < num_objects; num_object++) {
        radians[num_object] = 0;
    }

    setup_polygons(num_objects);

    draw_polygon(0, radians[0]);
    int last_num = 0;

    while (1) {
        int key = G_wait_key();
        int key_num = key - 49;

        if (key_num < num_objects) {
            G_rgb(1, 1, 1);
            G_fill_rectangle(0, 0, 800, 800);

            if (last_num == key_num) {
                radians[key_num] += (2 * M_PI) / 180;
                //printf("angle is now %lf\n", radians[key_num]);
            }
            last_num = key_num;

            draw_polygon(key_num, radians[key_num]);
        }
        else if (key == 'q') {
            G_save_image_to_file("t01c.xwd");
            G_close();
            exit(0) ;
        }
        else {
            continue;
        }
    }
    // hitting 0 should show the first image, and 1 the second, etc
}
