#include <bits/stdc++.h>

using namespace std;


#define cin(vec) for (auto &i : vec) cin >> i
#define cout(vec) for (auto &i : vec) cout << i << " ";cout << '\n';

#define pi (3.14159265358979323846)
#define EPS 0.01


bool AsmallerthanEPS (double a)
{
    return (a<EPS) && (-a < EPS);
}

vector<complex<double>> dft(vector<double>&signal)
{
    vector<complex<double>>ans;
    int N = signal.size();
    complex<double> i = complex<double>(0,1);
    // let zeta = e^(-2*pi*i/N);
    complex<double>zeta = exp(-2*pi/N*i);

    for(int k=0; k<N; k++)
    {
        complex<double>sum = 0;
        for(int n=0; n<N; n++)
        {
            sum += signal[n] * pow(zeta, n*k);
        }
        // rounding the real and the imaginary parts of the answer if they are so close to zero (compared using eps)
        sum = complex<double>((AsmallerthanEPS(real(sum))?0:real(sum)), (AsmallerthanEPS(imag(sum))?0:imag(sum)));
        
        ans.push_back(sum);
    }
    return ans;
}
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
    vector<double>signal(n);
    cin(signal);
    vector<complex<double>> x;
    for(auto& i:signal)
        x.push_back(complex<double>(i, 0));

    x = fft(x);
    cout(x);
}