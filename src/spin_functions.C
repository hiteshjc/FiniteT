#include"spin_functions.h"

double sz_fn(double spin,int ind)
{return double(ind)-spin;}

double sz2_fn(double spin,int ind)
{double sz=double(ind)-spin; return pow(sz,2.0);}

double splus_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin-sz)*(spin+sz+1)));}

double sminus_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin+sz)*(spin-sz+1)));}

double splus_sminus_fn(double spin,int ind)
{double sz=double(ind)-spin;return (spin+sz)*(spin-sz+1);}

double sminus_splus_fn(double spin,int ind)
{double sz=double(ind)-spin;return (spin-sz)*(spin+sz+1);}

double splus_splus_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin-sz)*(spin+sz+1)*(spin-sz-1)*(spin+sz+2)));}

double splus_sz_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin-sz)*(spin+sz+1)))*sz;}

double sz_splus_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin-sz)*(spin+sz+1)))*(sz+1);}

double sminus_sz_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin+sz)*(spin-sz+1)))*sz;}

double sz_sminus_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin+sz)*(spin-sz+1)))*(sz-1);}

double sminus_sminus_fn(double spin,int ind)
{double sz=double(ind)-spin;return sqrt(abs((spin+sz)*(spin-sz+1)*(spin+sz-1)*(spin-sz+2)));}
