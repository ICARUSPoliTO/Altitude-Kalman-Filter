#ifndef FILTER_C
#define FILTER_C

#include "filter.h"

//FILE *kf_log_fp = NULL;

void KalmanFilter_Init(KalmanFilter *kf)
{
    // Initialize state vector to zero
    zerosMatrix(kf->x, MAX_SIZE, MAX_SIZE);

    // Initial covariance: identity
    createIdentityMatrix(kf->P, MAX_SIZE);

    // Process noise matrix: diagonal
    float q_vec[MAX_SIZE] = {Q_POS, Q_VEL};
    createDiagonalMatrix(kf->Q, q_vec, MAX_SIZE);

    // Measurement matrix: [1 0]
    kf->H[0][0] = 1.0f;
    kf->H[0][1] = 0.0f;

    // Identity matrix
    createIdentityMatrix(kf->I, MAX_SIZE);
}

// check the alt data then create a csv file rahter than doing like this 
//          -> do file a operation and gather the data and do the kalman_filter thing on it then compare the csv files 
// pass alt data as parameter, add kf into (kf_altitude and kf_climbrate) the alt data struct 


/**
 * @brief Performs a single iteration of the Kalman filter update for altitude and vertical speed.
 * @param kf      Pointer to the KalmanFilter struct containing filter state and matrices.
 * @param z_meas  The latest measurement from the sensor.
 * @param dt      Time interval (in seconds) since the last update.
 */

void KalmanFilter_Update(KalmanFilter *kf, float z_meas,  float dt)
{
    float A[MAX_SIZE][MAX_SIZE] = {
        {1.0f, dt},
        {0.0f, 1.0f}
    };

    float x_pred[MAX_SIZE][1];
    float P_pred[MAX_SIZE][MAX_SIZE];

    // x_pred = A * x
    multiplyFloatMatrices(A, 2, 2, kf->x, 2, 1, x_pred);

    // P_pred = A * P * A' + Q
    float At[MAX_SIZE][MAX_SIZE];
    transposeFloatMatrix(A, At, 2, 2);

    float AP[MAX_SIZE][MAX_SIZE];
    multiplyFloatMatrices(A, 2, 2, kf->P, 2, 2, AP);

    float APA[MAX_SIZE][MAX_SIZE];
    multiplyFloatMatrices(AP, 2, 2, At, 2, 2, APA);

    addFloatMatrices(APA, kf->Q, P_pred, 2, 2);

    // Update step 
    // y = z_meas - H * x_pred
    float Hx[MAX_SIZE][1];
    multiplyFloatMatrices(kf->H, 1, 2, x_pred, 2, 1, Hx);
    float y = z_meas - Hx[0][0];

    // S = H * P_pred * H' + R
    float HP[MAX_SIZE][MAX_SIZE];
    multiplyFloatMatrices(kf->H, 1, 2, P_pred, 2, 2, HP);

    float HPHt[MAX_SIZE][MAX_SIZE];
    float Ht[MAX_SIZE][MAX_SIZE];
    transposeFloatMatrix(kf->H, Ht, 1, 2);
    multiplyFloatMatrices(HP, 1, 2, Ht, 2, 1, HPHt);

    float S = HPHt[0][0] + R;

    // K = P_pred * H' / S
    float PHt[MAX_SIZE][1];
    multiplyFloatMatrices(P_pred, 2, 2, Ht, 2, 1, PHt);

    float K[MAX_SIZE][1];
    for (int i = 0; i < 2; i++)
        K[i][0] = PHt[i][0] / S;

    // x = x_pred + K * y
    for (int i = 0; i < 2; i++)
        kf->x[i][1] = x_pred[i][0] + K[i][0] * y;

    // P = (I - K * H) * P_pred
    float KH[MAX_SIZE][MAX_SIZE];
    multiplyFloatMatrices(K, 2, 1, kf->H, 1, 2, KH);

    float IKH[MAX_SIZE][MAX_SIZE];
    subFloatMatrices(kf->I, KH, IKH, 2, 2);

    multiplyFloatMatrices(IKH, 2, 2, P_pred, 2, 2, kf->P);

}

#endif
