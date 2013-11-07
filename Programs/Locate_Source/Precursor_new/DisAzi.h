#ifndef DISAZI_H
#define DISAZI_H

template <class T>
int calc_dist (T lati1, T long1, T lati2, T long2, T *dist);

template <class T>
int calc_azimuth(T lati1, T long1, T lati2, T long2, T *alpha1);

#endif
