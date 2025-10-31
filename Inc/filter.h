#ifndef FILTER_H
#define FILTER_H

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

#define R 0.4f // measurement noise variance (m^2) defualt 0.4
#define Q_POS 1e-5f // process noise (altitude) default 1e-5
#define Q_VEL 1e-3f // process noise (velocity) default 1e-3

typedef struct {
    float x[MAX_SIZE][1];                   // States matrix: x = [altitude; vertical_velocity]
    float P[MAX_SIZE][MAX_SIZE];         // (Initial) Covariance matrix
    float Q[MAX_SIZE][MAX_SIZE];         // Process noise matrix
    float H[1][MAX_SIZE];                // Measurement matrix
    float I[MAX_SIZE][MAX_SIZE];         // Identity matrix
} KalmanFilter;

void KalmanFilter_Init(KalmanFilter *kf);
void KalmanFilter_Update(KalmanFilter *kf, float z_meas, float dt);

#endif 