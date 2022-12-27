#include <bits/stdc++.h>
#include "AudioFile.h"

#define pi (3.14159265358979323846)
#define cin(vec) for (auto &i : vec) cin >> i
#define cout(vec) for (auto &i : vec) cout << i << " ";cout << '\n';

using namespace std;
using cd = complex<double>;

const double PI = acos(-1);


vector<complex<double>> fft(vector<complex<double>>& p)
{
    int n = p.size();
    if(n==1){
        return p;
    }
    int nn = pow(2, (int)(ceil(log2(p.size()))));   // smallest power of two smaller than or equal to length of p 
    complex<double>i = complex<double>(0, 1);
    complex<double>omega = exp(2* pi * i / (double)n); // Idft omega = omega * (-1/n)
    vector<complex<double>>pe(n/2), po(n/2);
    for(int j=0; j<n; j++)
    {
        if(j%2==0)
            pe[j/2] = p[j];
        else
            po[j/2] = p[j];
    }
    vector<complex<double>> ye = fft(pe);
    vector<complex<double>> yo = fft(po);
    vector<complex<double>> y(n);
    double ang = 2 * pi / n;
    complex<double> w(1), wn(cos(ang), sin(ang));
    for(int j=0; j<n/2; j++)
    {
        y[j] = ye[j] + yo[j]*w;
        y[j+n/2] = ye[j] - yo[j]*w;
        w *= wn;
    }
    return y;
}

int main()
{
    int n;      cin>>n;
    vector<complex<double>>signal(n);
    for(int i=0; i<n; i++)
    {
        int x;      cin>>x;
        signal[i] = complex<double>(x, 0);
    }
    vector<complex<double>>ans = fft(signal);
    for(auto& i:ans)
    {
        cout << i << " ";
    }
}