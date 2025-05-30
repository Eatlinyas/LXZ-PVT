/** Modified Dip Latitude grid.
 *
 * @author Angela Aragon-Angel (maria-angeles.aragon@ec.europa.eu)
 * @ingroup NeQuickG_JRC
 * @copyright Joint Research Centre (JRC), 2019<br>
 *  This software has been released as free and open source software
 *  under the terms of the European Union Public Licence (EUPL), version 1.<br>
 *  Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
 * @file
 */
#include "NeQuickG_JRC_MODIP_grid.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "NeQuickG_JRC.h"
#ifndef FTR_MODIP_CCIR_AS_CONSTANTS
#include "NeQuickG_JRC_error.h"
#endif //FTR_MODIP_CCIR_AS_CONSTANTS
#include "NeQuickG_JRC_interpolate.h"
#include "NeQuickG_JRC_macros.h"

#define NEQUICK_G_JRC_MODIP_GRID_LONG_UNIQUE_COUNT \
  (360/NEQUICK_G_JRC_MODIP_GRID_LONG_STEP_DEGREE)

/** See F.2.4.1.
 *  In the pseudocode this correction offset is present
 *  But ESA, in its public version, do not use it
 *  and adds an explanation in the code.
 *  we follow the pseudo-code.
 */
#define NEQUICK_G_JRC_MODIP_GRID_LAT_OFFSET_CORRECTION (1.0e-6)

#ifdef FTR_MODIP_CCIR_AS_CONSTANTS
/** MODIP matrix available at its first release dated on year 2001 and calculated at a
 * height of 300 km (Private communication, Galileo Service Centre, 2018).
 * Private communication, Galileo Service Centre, ID ticket #569, 2018 */
static const double_t g_corrected_modip_degree
    [NEQUICK_G_JRC_MODIP_GRID_LAT_POINTS_COUNT]
    [NEQUICK_G_JRC_MODIP_GRID_LONG_POINTS_COUNT] = {
  { -76.37, -76.34, -76.30, -76.25, -76.55, -76.47, -76.75, -76.65, -76.90, -76.79, -77.04, -77.27, -77.19, -77.44, -77.37, -77.33, -77.30, -77.30, -77.31, -77.34, -77.40, -77.13, -77.21, -76.96, -77.06, -76.81, -76.91, -76.66, -76.76, -76.48, -76.56, -76.25, -76.31, -76.34, -76.37, -76.37, -76.37, -76.34, -76.30 },
  { -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00, -90.00 },
  { -77.31, -77.34, -77.40, -77.13, -77.21, -76.96, -77.06, -76.81, -76.91, -76.66, -76.76, -76.48, -76.56, -76.25, -76.31, -76.34, -76.37, -76.37, -76.37, -76.34, -76.30, -76.25, -76.55, -76.47, -76.75, -76.65, -76.90, -76.79, -77.04, -77.27, -77.19, -77.44, -77.37, -77.33, -77.30, -77.30, -77.31, -77.34, -77.40 },
  { -72.98, -73.10, -72.84, -72.61, -72.40, -72.20, -72.00, -71.80, -71.59, -71.36, -71.11, -70.82, -70.49, -70.64, -70.23, -70.31, -70.35, -70.36, -70.34, -70.28, -70.18, -70.59, -70.42, -70.74, -71.02, -71.26, -71.47, -71.66, -72.30, -72.47, -72.65, -72.86, -73.10, -73.38, -73.32, -73.32, -72.98, -73.10, -72.84 },
  { -70.32, -70.13, -69.97, -69.37, -69.25, -68.64, -67.99, -67.89, -67.19, -67.04, -66.22, -65.96, -65.61, -65.84, -65.33, -65.41, -65.43, -65.42, -65.35, -65.26, -65.78, -65.58, -65.97, -66.28, -66.51, -67.28, -67.95, -68.01, -68.58, -69.12, -69.65, -70.20, -70.34, -70.58, -70.51, -70.58, -70.32, -70.13, -69.97 },
  { -67.88, -67.33, -66.76, -66.71, -66.11, -65.50, -64.87, -64.23, -63.55, -62.82, -62.72, -61.76, -61.41, -60.88, -61.01, -61.01, -60.95, -60.84, -61.51, -61.36, -61.95, -61.71, -62.13, -62.42, -63.30, -64.03, -64.63, -65.15, -66.19, -66.60, -67.52, -67.93, -68.39, -68.55, -68.58, -68.42, -67.88, -67.33, -66.76 },
  { -65.43, -64.90, -64.33, -63.71, -63.06, -62.40, -61.75, -61.11, -60.47, -59.79, -58.18, -58.04, -57.67, -57.02, -57.01, -57.72, -57.46, -58.06, -57.84, -58.51, -58.32, -58.92, -59.38, -59.66, -60.54, -61.22, -61.75, -62.83, -63.80, -64.68, -65.51, -66.31, -66.61, -67.17, -66.81, -65.90, -65.43, -64.90, -64.33 },
  { -63.10, -62.52, -61.89, -61.21, -60.49, -59.77, -59.08, -58.42, -57.00, -56.36, -54.71, -53.70, -53.29, -53.45, -53.20, -53.67, -55.00, -55.41, -55.94, -55.71, -55.57, -56.29, -55.98, -57.17, -58.06, -58.70, -59.90, -60.92, -61.83, -62.70, -63.55, -64.45, -64.56, -64.70, -64.14, -63.63, -63.10, -62.52, -61.89 },
  { -60.55, -59.89, -59.18, -58.41, -57.60, -56.79, -56.00, -55.29, -53.77, -53.17, -51.50, -50.51, -50.08, -50.08, -50.56, -51.68, -52.65, -52.78, -54.07, -53.77, -53.68, -54.54, -54.38, -54.82, -55.78, -56.42, -57.62, -58.65, -60.26, -61.19, -61.61, -62.29, -62.19, -61.90, -61.80, -61.18, -60.55, -59.89, -59.18 },
  { -57.52, -56.76, -56.72, -55.88, -54.99, -54.10, -53.23, -52.44, -50.80, -49.16, -48.52, -46.43, -45.95, -45.73, -47.00, -48.82, -50.45, -51.19, -52.26, -52.75, -52.68, -52.80, -52.84, -53.46, -53.66, -54.34, -55.58, -57.41, -58.42, -58.87, -59.57, -59.97, -60.07, -59.79, -58.99, -58.25, -57.52, -56.76, -56.72 },
  { -55.35, -54.54, -53.68, -52.76, -51.79, -50.80, -49.83, -48.94, -48.17, -46.46, -44.70, -42.50, -41.96, -42.75, -44.84, -46.20, -48.49, -49.87, -51.57, -51.85, -51.71, -51.95, -51.34, -51.31, -52.58, -53.34, -54.63, -55.78, -56.24, -57.62, -57.94, -57.92, -57.47, -57.17, -56.26, -56.16, -55.35, -54.54, -53.68 },
  { -52.61, -51.73, -50.81, -49.82, -48.79, -47.75, -46.73, -45.73, -44.88, -43.03, -41.20, -38.94, -38.31, -38.86, -41.88, -43.95, -46.94, -48.96, -50.32, -51.24, -51.85, -51.24, -50.87, -50.24, -50.85, -51.80, -52.39, -53.74, -54.41, -55.38, -55.25, -55.40, -54.97, -54.65, -53.63, -53.49, -52.61, -51.73, -50.81 },
  { -50.08, -48.20, -47.20, -47.17, -46.12, -45.07, -42.87, -41.78, -40.82, -38.82, -36.89, -34.51, -33.78, -35.48, -38.12, -41.01, -44.83, -47.57, -49.58, -50.20, -50.62, -50.87, -49.77, -49.45, -49.44, -49.72, -51.40, -52.11, -53.01, -52.69, -52.77, -53.07, -52.64, -51.44, -51.17, -50.13, -50.08, -48.20, -47.20 },
  { -46.90, -44.89, -43.84, -42.71, -41.57, -40.48, -39.39, -36.91, -35.88, -33.72, -31.73, -29.23, -28.33, -29.76, -33.48, -38.67, -42.12, -45.72, -47.53, -48.89, -50.02, -50.20, -49.21, -48.22, -47.56, -48.11, -49.16, -50.11, -50.41, -50.29, -50.57, -50.12, -49.68, -49.27, -48.03, -47.90, -46.90, -44.89, -43.84 },
  { -42.94, -40.78, -39.69, -38.48, -37.26, -36.16, -33.63, -32.45, -29.87, -27.64, -23.91, -21.22, -21.91, -22.92, -29.50, -34.27, -40.00, -43.37, -46.09, -47.36, -48.44, -48.57, -47.60, -45.76, -45.35, -46.19, -46.51, -46.76, -47.29, -47.35, -46.83, -46.42, -45.93, -45.45, -45.13, -43.97, -42.94, -40.78, -39.69 },
  { -38.00, -37.05, -34.60, -33.29, -32.00, -29.31, -28.16, -25.31, -22.60, -20.39, -16.58, -13.78, -12.32, -16.69, -22.92, -30.45, -35.94, -40.47, -43.26, -45.67, -45.82, -45.96, -44.97, -43.15, -42.92, -42.92, -43.46, -43.88, -44.62, -43.76, -43.30, -42.89, -42.33, -41.78, -40.28, -39.04, -38.00, -37.05, -34.60 },
  { -33.30, -30.99, -28.38, -26.97, -25.53, -22.62, -19.66, -18.45, -15.74, -11.72, -7.87, -4.92, -3.12, -7.12, -14.93, -24.08, -32.47, -36.99, -40.00, -41.64, -43.06, -42.18, -41.13, -40.42, -39.12, -38.01, -38.74, -39.34, -40.26, -39.38, -38.94, -38.49, -37.82, -37.16, -35.51, -34.24, -33.30, -30.99, -28.38 },
  { -25.81, -23.45, -20.73, -19.23, -15.79, -14.45, -11.27, -8.07, -5.36, -1.38, 1.54, 4.46, 4.34, -0.15, -7.30, -18.17, -26.70, -32.91, -36.30, -38.36, -38.91, -38.05, -36.96, -34.95, -33.56, -32.35, -33.22, -33.93, -35.02, -34.07, -33.58, -33.04, -32.18, -31.36, -29.52, -28.25, -25.81, -23.45, -20.73 },
  { -18.21, -14.07, -11.29, -9.72, -7.95, -4.47, -1.12, -0.13, 3.20, 7.24, 11.16, 13.63, 13.20, 8.68, -0.02, -8.99, -19.80, -26.62, -30.66, -31.94, -32.88, -32.04, -29.36, -27.06, -25.49, -25.75, -26.72, -27.51, -27.17, -27.63, -25.41, -24.66, -25.14, -24.09, -22.00, -20.69, -18.21, -14.07, -11.29 },
  { -8.65, -4.46, -1.73, 0.13, 2.22, 4.54, 7.05, 9.67, 12.82, 16.46, 19.79, 21.67, 21.01, 16.82, 8.72, -0.07, -11.56, -19.23, -22.42, -24.32, -24.10, -23.26, -20.19, -17.56, -15.77, -16.03, -17.09, -17.97, -17.55, -17.99, -17.17, -16.10, -16.39, -14.98, -14.45, -11.16, -8.65, -4.46, -1.73 },
  { 0.78, 4.62, 7.90, 10.21, 12.09, 14.13, 16.34, 18.68, 21.42, 24.45, 27.06, 28.37, 27.58, 23.89, 16.91, 7.39, -1.72, -8.45, -12.47, -13.20, -13.31, -12.51, -9.01, -7.96, -5.94, -4.22, -5.36, -6.30, -7.78, -6.17, -5.05, -5.58, -5.56, -5.74, -4.93, -1.35, 0.78, 4.62, 7.90 },
  { 10.53, 14.10, 17.12, 19.21, 20.83, 22.54, 24.39, 26.39, 28.66, 31.06, 33.02, 33.86, 32.99, 29.88, 24.18, 16.35, 8.11, 1.73, -0.45, -1.65, -0.05, -0.74, 1.49, 4.17, 5.96, 6.21, 5.41, 4.44, 3.94, 4.19, 4.79, 4.99, 4.61, 4.42, 5.26, 7.37, 10.53, 14.10, 17.12 },
  { 19.28, 22.29, 24.88, 26.70, 28.07, 29.49, 31.06, 32.75, 34.60, 36.47, 37.90, 38.39, 37.50, 34.90, 30.42, 24.42, 18.03, 12.98, 10.38, 9.92, 10.50, 11.69, 13.59, 15.77, 17.22, 17.43, 16.81, 16.02, 15.59, 15.74, 16.07, 15.96, 15.30, 14.76, 15.12, 16.69, 19.28, 22.29, 24.88 },
  { 26.70, 29.05, 31.16, 32.70, 33.90, 35.12, 36.47, 37.92, 39.45, 40.91, 41.96, 42.19, 41.32, 39.16, 35.70, 31.30, 26.70, 23.08, 21.25, 21.01, 21.61, 22.63, 24.09, 25.69, 26.74, 26.93, 26.51, 25.95, 25.66, 25.72, 25.83, 25.49, 24.69, 23.90, 23.83, 24.79, 26.70, 29.05, 31.16 },
  { 32.85, 34.58, 36.23, 37.53, 38.62, 39.71, 40.91, 42.19, 43.48, 44.66, 45.43, 45.49, 44.65, 42.86, 40.21, 37.04, 33.89, 31.44, 30.23, 30.14, 30.66, 31.45, 32.51, 33.62, 34.37, 34.53, 34.28, 33.93, 33.76, 33.79, 33.75, 33.31, 32.45, 31.54, 31.14, 31.57, 32.85, 34.58, 36.23 },
  { 37.97, 39.18, 40.43, 41.54, 42.53, 43.56, 44.66, 45.80, 46.94, 47.90, 48.48, 48.44, 47.65, 46.16, 44.13, 41.87, 39.72, 38.12, 37.34, 37.33, 37.74, 38.35, 39.10, 39.85, 40.37, 40.53, 40.41, 40.25, 40.16, 40.17, 40.07, 39.59, 38.73, 37.78, 37.19, 37.23, 37.97, 39.18, 40.43 },
  { 42.38, 43.15, 44.10, 45.03, 45.95, 46.92, 47.95, 49.00, 50.00, 50.82, 51.27, 51.16, 50.44, 49.18, 47.61, 45.97, 44.51, 43.46, 42.97, 42.99, 43.30, 43.76, 44.30, 44.81, 45.18, 45.33, 45.33, 45.29, 45.30, 45.31, 45.15, 44.67, 43.86, 42.93, 42.23, 42.02, 42.38, 43.15, 44.10 },
  { 46.29, 46.76, 47.46, 48.25, 49.09, 50.00, 50.97, 51.93, 52.82, 53.52, 53.87, 53.73, 53.07, 52.02, 50.78, 49.57, 48.55, 47.84, 47.52, 47.55, 47.79, 48.13, 48.51, 48.87, 49.16, 49.32, 49.40, 49.46, 49.53, 49.55, 49.38, 48.91, 48.15, 47.29, 46.56, 46.20, 46.29, 46.76, 47.46 },
  { 49.92, 50.19, 50.69, 51.34, 52.10, 52.95, 53.84, 54.71, 55.50, 56.10, 56.38, 56.21, 55.61, 54.73, 53.75, 52.83, 52.07, 51.57, 51.35, 51.37, 51.54, 51.80, 52.07, 52.34, 52.58, 52.76, 52.91, 53.05, 53.17, 53.19, 53.02, 52.58, 51.89, 51.10, 50.41, 49.98, 49.92, 50.19, 50.69 },
  { 53.37, 53.50, 53.87, 54.42, 55.09, 55.84, 56.64, 57.42, 58.11, 58.62, 58.82, 58.64, 58.11, 57.37, 56.57, 55.85, 55.28, 54.90, 54.72, 54.72, 54.84, 55.02, 55.22, 55.44, 55.65, 55.86, 56.06, 56.26, 56.41, 56.45, 56.28, 55.87, 55.26, 54.56, 53.93, 53.52, 53.37, 53.50, 53.87 },
  { 56.74, 56.80, 57.06, 57.50, 58.07, 58.72, 59.42, 60.10, 60.69, 61.11, 61.25, 61.05, 60.58, 59.96, 59.31, 58.73, 58.28, 57.98, 57.82, 57.79, 57.86, 57.98, 58.14, 58.33, 58.53, 58.75, 59.00, 59.23, 59.41, 59.46, 59.31, 58.95, 58.41, 57.82, 57.28, 56.90, 56.74, 56.80, 57.06 },
  { 60.07, 60.09, 60.28, 60.62, 61.09, 61.63, 62.21, 62.78, 63.28, 63.61, 63.69, 63.48, 63.06, 62.54, 62.02, 61.56, 61.19, 60.93, 60.79, 60.74, 60.77, 60.85, 60.98, 61.14, 61.34, 61.58, 61.83, 62.09, 62.28, 62.34, 62.21, 61.91, 61.47, 60.98, 60.54, 60.22, 60.07, 60.09, 60.28 },
  { 63.41, 63.42, 63.55, 63.81, 64.16, 64.59, 65.05, 65.50, 65.91, 66.17, 66.19, 65.97, 65.60, 65.17, 64.76, 64.39, 64.10, 63.88, 63.75, 63.69, 63.70, 63.75, 63.85, 64.00, 64.19, 64.42, 64.67, 64.92, 65.11, 65.17, 65.08, 64.84, 64.50, 64.13, 63.79, 63.54, 63.41, 63.42, 63.55 },
  { 66.83, 66.83, 66.92, 67.10, 67.35, 67.65, 68.00, 68.34, 68.64, 68.83, 68.80, 68.57, 68.26, 67.93, 67.61, 67.34, 67.11, 66.94, 66.82, 66.77, 66.76, 66.80, 66.88, 67.01, 67.18, 67.38, 67.60, 67.82, 67.99, 68.07, 68.02, 67.85, 67.60, 67.34, 67.10, 66.93, 66.83, 66.83, 66.92 },
  { 70.42, 70.42, 70.48, 70.58, 70.74, 70.94, 71.16, 71.40, 71.61, 71.72, 71.63, 71.42, 71.18, 70.94, 70.72, 70.53, 70.37, 70.25, 70.16, 70.12, 70.11, 70.14, 70.20, 70.30, 70.43, 70.58, 70.76, 70.93, 71.07, 71.15, 71.14, 71.04, 70.90, 70.74, 70.59, 70.48, 70.42, 70.42, 70.48 },
  { 74.37, 74.37, 74.40, 74.46, 74.54, 74.65, 74.77, 74.90, 75.03, 75.03, 74.90, 74.75, 74.60, 74.46, 74.33, 74.21, 74.12, 74.05, 73.99, 73.97, 73.96, 73.98, 74.02, 74.08, 74.16, 74.26, 74.36, 74.47, 74.57, 74.64, 74.66, 74.64, 74.58, 74.52, 74.45, 74.40, 74.37, 74.37, 74.40 },
  { 79.15, 79.15, 79.17, 79.19, 79.22, 79.25, 79.29, 79.31, 79.28, 79.23, 79.18, 79.12, 79.06, 79.01, 78.96, 78.92, 78.88, 78.85, 78.83, 78.82, 78.82, 78.83, 78.84, 78.87, 78.90, 78.93, 78.97, 79.02, 79.06, 79.10, 79.13, 79.15, 79.16, 79.16, 79.15, 79.15, 79.15, 79.15, 79.17 },
  { 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00, 90.00 },
  { 78.83, 78.82, 78.82, 78.83, 78.84, 78.87, 78.90, 78.93, 78.97, 79.02, 79.06, 79.10, 79.13, 79.15, 79.16, 79.16, 79.15, 79.15, 79.15, 79.15, 79.17, 79.19, 79.22, 79.25, 79.29, 79.31, 79.28, 79.23, 79.18, 79.12, 79.06, 79.01, 78.96, 78.92, 78.88, 78.85, 78.83, 78.82, 78.82 },
};
#endif // FTR_MODIP_CCIR_AS_CONSTANTS

typedef struct grid_position_st {
  uint8_t index;
  double_t offset;
} grid_position_t;

#ifndef FTR_MODIP_CCIR_AS_CONSTANTS
static int32_t grid_file_open(
  const char* const pFileName,
  FILE** ppFile) {

  *ppFile = NULL;
  if (!pFileName) {
    NEQUICK_ERROR_RETURN(
      NEQUICK_ERROR_SRC_MODIP,
      NEQUICK_ERROR_CODE_BAD_MODIP_GRID_FILE,
      "Modip file. null pointer exception.");
  }
  // to satisfy lint
  assert(pFileName);
  *ppFile = fopen(pFileName, "r");
  if (!(*ppFile)) {
    NEQUICK_ERROR_RETURN(
      NEQUICK_ERROR_SRC_MODIP,
      NEQUICK_ERROR_CODE_BAD_MODIP_GRID_FILE,
      "Cannot open Modip file %s", pFileName);
  }
  return NEQUICK_OK;
}

static int32_t grid_load(
  FILE * const pGridFile,
  modip_grid_t* const pGrid) {

  size_t latitude_index;
  size_t longitude_index;
  for (latitude_index = 0;
     latitude_index < NEQUICK_G_JRC_MODIP_GRID_LAT_POINTS_COUNT;
     latitude_index++) {
    for (longitude_index = 0;
       longitude_index < NEQUICK_G_JRC_MODIP_GRID_LONG_POINTS_COUNT;
       longitude_index++) {
      if (fscanf(pGridFile, "%lf",
        &pGrid->corrected_modip_degree
          [latitude_index][longitude_index])
          != 0x01) {
        NEQUICK_ERROR_RETURN(
          NEQUICK_ERROR_SRC_MODIP,
          NEQUICK_ERROR_CODE_BAD_MODIP_STRUCTURE,
          "error scanning modip grid using fscanf");
      }
    }
  }
  pGrid->is_loaded = true;
  return NEQUICK_OK;
}

int32_t modip_grid_init(
  const char* const pModip_file,
  modip_grid_t* const pGrid) {

  int32_t ret;

  pGrid->is_loaded = false;
  FILE* pGrid_file;
  ret = grid_file_open(pModip_file, &pGrid_file);
  if (ret != NEQUICK_OK) {
    return ret;
  }

  ret = grid_load(pGrid_file, pGrid);

  if (pGrid_file) {
    fclose(pGrid_file);
  }

  return ret;
}
#endif // FTR_MODIP_CCIR_AS_CONSTANTS

/** TBD: in the ESA version this uses a while
 * so any value is finally corrected.
 */
static void grid_long_fix_index(
  grid_position_t * const pGrid_position, int8_t index) {
  if (index < 0) {
    assert(index > NEQUICK_G_JRC_MODIP_GRID_LONG_UNIQUE_COUNT);
    pGrid_position->index = (uint8_t)(index +
      NEQUICK_G_JRC_MODIP_GRID_LONG_UNIQUE_COUNT);
  } else if (index >=
      NEQUICK_G_JRC_MODIP_GRID_LONG_UNIQUE_COUNT) {
    pGrid_position->index = (uint8_t)(index -
      NEQUICK_G_JRC_MODIP_GRID_LONG_UNIQUE_COUNT);
  } else {
    pGrid_position->index = (uint8_t)index;
  }
}

static void grid_get_long_position(
  double_t longitude_degree,
  grid_position_t * const pGrid_position) {

  double_t index_with_offset =
    (longitude_degree + NEQUICK_G_JRC_LONG_MAX_VALUE_DEGREE) /
    NEQUICK_G_JRC_MODIP_GRID_LONG_STEP_DEGREE;

  double_t _index = floor(index_with_offset);
  pGrid_position->offset = (index_with_offset - _index);

  // and, if required, adjust for sign and wrap to grid
  grid_long_fix_index(pGrid_position, (int8_t)_index);
}

static void grid_get_lat_position(
  double_t latitude_degree,
  grid_position_t * const pGrid_position) {

  double_t index_with_offset =
    (latitude_degree + NEQUICK_G_JRC_LAT_MAX_VALUE_DEGREE) /
    NEQUICK_G_JRC_MODIP_GRID_LAT_STEP_DEGREE;

  assert(index_with_offset > 0.0);
  double_t _index =
    floor(index_with_offset -
                   NEQUICK_G_JRC_MODIP_GRID_LAT_OFFSET_CORRECTION);
  pGrid_position->index = (uint8_t)_index;
  pGrid_position->offset = (index_with_offset - pGrid_position->index);
}

#ifdef FTR_MODIP_CCIR_AS_CONSTANTS
#define grid_get_lon_interpol_point(pGrid, pLatitude, longitude_grid_index) \
  grid_get_lon_interpol_point(pLatitude, longitude_grid_index)
#endif

/**
 * Interpolate across lat grid
 */
static double_t
grid_get_lon_interpol_point(
  const modip_grid_t* const pGrid,
  const grid_position_t* const pLatitude,
  uint_fast8_t longitude_grid_index) {
  double_t
    lat_interpol_points[NEQUICK_G_JRC_INTERPOLATE_POINT_COUNT];
  size_t lat_grid_index = pLatitude->index;
  for (size_t i = 0;
        i < NEQUICK_G_JRC_INTERPOLATE_POINT_COUNT; i++) {
    lat_interpol_points[i] =
#ifdef FTR_MODIP_CCIR_AS_CONSTANTS
      g_corrected_modip_degree
#else
      pGrid->corrected_modip_degree
#endif
        [lat_grid_index++]
        [longitude_grid_index];
  }
  return interpolation_third_order(
    lat_interpol_points,
    pLatitude->offset);
}

#ifdef FTR_MODIP_CCIR_AS_CONSTANTS
#define grid_2D_interpolation(pGrid, pLongitude, pLatitude) \
  grid_2D_interpolation(pLongitude, pLatitude)
#endif

static double_t grid_2D_interpolation(
  const modip_grid_t* const pGrid,
  grid_position_t* pLongitude,
  const grid_position_t* const pLatitude) {

  double_t
    lon_interpol_points[NEQUICK_G_JRC_INTERPOLATE_POINT_COUNT];

  for (size_t i = 0;
        i < NEQUICK_G_JRC_INTERPOLATE_POINT_COUNT; i++) {
    lon_interpol_points[i] = grid_get_lon_interpol_point(
                              pGrid, pLatitude, pLongitude->index++);
  }
  return interpolation_third_order(
    lon_interpol_points,
    pLongitude->offset);
}

void modip_grid_interpolate(
  const modip_grid_t* const pGrid,
  const position_t* const pPosition,
  double_t *pModip_degree) {

  *pModip_degree = 0;

  // Obtain grid <index, offset> for longitude
  grid_position_t longitude;
  grid_get_long_position(pPosition->longitude.degree, &longitude);

  // Obtain grid <index, offset> for latitude
  grid_position_t latitude;
  grid_get_lat_position(pPosition->latitude.degree, &latitude);

  *pModip_degree = grid_2D_interpolation(pGrid, &longitude, &latitude);
}

#undef NEQUICK_G_JRC_MODIP_GRID_LONG_UNIQUE_COUNT
