#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <FPT.h>
#include "M2d_matrix_toolsS.c"
#include "intersect_two.c"

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


double clicked_points[2][100];
int clicked_points_count;
double clicked_center[2];

void find_clicked_center() {
    clicked_center[0] = 0;
    clicked_center[1] = 0;

    for (int point = 0; point < clicked_points_count; point++) {
        clicked_center[0] += clicked_points[0][point];
        clicked_center[1] += clicked_points[1][point];
    }

    clicked_center[0] /= clicked_points_count;
    clicked_center[1] /= clicked_points_count;
}

void click_polygon() {
    clicked_points_count = 0;

    G_rgb(0, 1, 0);
    while (1) {
        double p[2];
        G_wait_click(p);

        if (p[0] < 50 && p[1] < 50) {
            break;
        }

        G_circle(p[0], p[1], 5);

        clicked_points[0][clicked_points_count] = p[0];
        clicked_points[1][clicked_points_count] = p[1];
        clicked_points_count++;
    }

    find_clicked_center();
}

int in_cut(double p1[2], double p2[2], double p3[2]) {
    double x_factor, y_factor, intercept;

    // y = mx + b
    // 0 = mx - 1y + b
    // 0 = ax + By + C
    // vertical   lines have x = k, 1x + 0y - k = 0
    // horizontal lines have y = k, 0x + 1y - k = 0
    if (p2[0] == p1[0]) { // if vertical
        x_factor = 1;
        y_factor = 0;
        intercept = -1 * p1[0];
    }

    else if (p2[1] == p1[1]) { // if horizontal
        x_factor = 0;
        y_factor = 1;
        intercept = -1 * p1[1];
    }
    else { // if normal line
        x_factor = (p2[1] - p1[1]) / (p2[0] - p1[0]);
        y_factor = -1;
        intercept = p1[1] - (p1[0] * x_factor);
    }

    // equation Ax + By + C = 0
    double center_sign = (x_factor * clicked_center[0]) + (y_factor * clicked_center[1]) + intercept;
    double point_sign = (x_factor * p3[0]) + (y_factor * p3[1]) + intercept;

    if (center_sign > 0 && point_sign > 0) {
        return 0; // on same side as center
    }
    else if (center_sign < 0 && point_sign < 0) {
        return 0; // on same side as center
    }

    return 1; // on other side of center
}



// do iterative infinite slices off of the figure
// for each line in the clicked_polygon, cut an infinite line off of the figure
// check for each point, is it below/above the line?
// if so, no change
// if not, find intersection
// will probably give you more points than you started with!
void cut_polygon(double poly_points[2][100], int length) {
    double c1[2], c2[2], new_points[2][100], p1[2], p2[2], intersections[2];
    int check_intersection, p1_cut, p2_cut;
    int counter = 0;

    for (int clicked_line = 0; clicked_line < clicked_points_count; clicked_line++) {
        // sets the x and y values for v1

        c1[0] = clicked_points[0][clicked_line];
        c1[1] = clicked_points[1][clicked_line];

        // sets c2 values to next point, unless vertex is the last in array

        if (clicked_line < (clicked_points_count - 1)) {
            c2[0] = clicked_points[0][clicked_line + 1];
            c2[1] = clicked_points[1][clicked_line + 1];
        }
        else {
            c2[0] = clicked_points[0][0];
            c2[1] = clicked_points[1][0];
        }

	    counter = 0;

        for (int point = 0; point < length; point++) {
            // are the x values above the center and below the line?

            p1[0] = poly_points[0][point];
            p1[1] = poly_points[1][point];

            if (point < (length - 1)) {
                p2[0] = poly_points[0][point + 1];
                p2[1] = poly_points[1][point + 1];
            }
            else {
                p2[0] = poly_points[0][0];
                p2[1] = poly_points[1][0];
            }

            // check if point is in the cut
    	    p1_cut = in_cut(c1, c2, p1);
    	    p2_cut = in_cut(c1, c2, p2);


    	    // if both good
    	    if (p1_cut == 0 && p2_cut == 0) {
                new_points[0][counter] = p2[0];
                new_points[1][counter] = p2[1];
                counter++;
    	    }
    	    // if just one good
    	    // only need to check p1, as p2 will be found next
    	    else if (p1_cut == 1 && p2_cut == 0) {

                check_intersection = intersect_2_lines(c1, c2, p1, p2, intersections);
                    new_points[0][counter] = intersections[0];
                    new_points[1][counter] = intersections[1];
                    new_points[0][counter + 1] = p2[0];
                    new_points[1][counter + 1] = p2[1];
                    counter++;
                    counter++;
            }
            // if first good
            else if (p1_cut == 0 && p2_cut == 1) {

                check_intersection = intersect_2_lines(c1, c2, p1, p2, intersections);
                    new_points[0][counter] = intersections[0];
                    new_points[1][counter] = intersections[1];
                    counter++;
            }


        } // for point

        // set poly_points to new_points
    	// set length to counter
    	for (int copy = 0; copy < counter; copy++) {
      	    poly_points[0][copy] = new_points[0][copy];
    	    poly_points[1][copy] = new_points[1][copy];
    	}

    	length = counter;
        counter = 0;
    } // for clicked

    G_fill_polygon(poly_points[0], poly_points[1], length);
}



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

void draw_polygon(int file_number, double radians, int number) {
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

        if (number == 1) {
            cut_polygon(poly_points, length);
        }
        else if (number == 0) {
	        G_fill_polygon(poly_points[0], poly_points[1], length);
        }
    }
}

void setup_polygons(int num_objects) {
    for (int file_number = 0; file_number < num_objects; file_number++) {
        find_center(file_number);
        find_scale(file_number);
    }
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

        // next NUMPOINTS lines are the x and y values of each point
		for (int i = 0; i < num_points[k]; i++) {
			fscanf(fp, "%lf %lf", &x[k][i], &y[k][i]);
		}

        // next line is the number of polygons
		fscanf(fp, "%d", &num_polys[k]);

        // next NUMPOLY lines are the number of lines in this polygon
		for (int i = 0; i < num_polys[k]; i++) {
			fscanf(fp, "%d", &psize[k][i]);

            // followed by the points that are connected in this polygon
			for (int j = 0; j < psize[k][i]; j++) {
				fscanf(fp, "%d", &con[k][i][j]);
			}
		}

        // next NUMPOLY lines are the RGB values of each polygon
		for (int i = 0; i < num_polys[k]; i++) {
			fscanf(fp, "%lf %lf %lf", &red[k][i], &grn[k][i], &blu[k][i]);
		}
	}

	// draw these polygons!
    G_init_graphics(600, 600);
    double radians[num_objects];

    for (int num_object = 0; num_object < num_objects; num_object++) {
        radians[num_object] = 0;
    }

    setup_polygons(num_objects);


    G_rgb(0, 0, 0);
    G_fill_rectangle(0, 0, 600, 600);
    G_rgb(0, 1, 0);
    G_fill_rectangle(0, 0, 50, 50);
    draw_polygon(0, radians[0], 0);

    int last_num = 0;

    click_polygon();
    G_rgb(0, 0, 0);
    G_fill_rectangle(0, 0, 600, 600);

    draw_polygon(0, radians[0], 1);

    G_rgb(1, 1, 1);
    G_fill_polygon(clicked_points[0], clicked_points[1], clicked_points_count);
    draw_polygon(0, radians[0], 1); // draw again as first draw will be overwritten

    while (1) {
        G_rgb(0, 0, 1);
        G_polygon(clicked_points[0], clicked_points[1], clicked_points_count);
        int key = G_wait_key();
        int key_num = key - 49;

        if (key_num < num_objects) {
            G_rgb(0, 0, 0);
            G_fill_rectangle(0, 0, 600, 600);

            G_rgb(1, 1, 1);
            G_fill_polygon(clicked_points[0], clicked_points[1], clicked_points_count);

            if (last_num == key_num) {
                radians[key_num] += (2 * M_PI) / 180;
            }
            last_num = key_num;

            draw_polygon(key_num, radians[key_num], 1);
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
