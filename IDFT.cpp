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

complex<double>multiplycomplex(complex<double>a, complex<double>b)
{
    double r1,i1,r2,i2;
    r1 = real(a);
    i1 = imag(a);
    r2 = real(b);
    i2 = imag(b);
    return complex<double>((r1*r2-i1*i2), (r1*i2 + r2*i1));
}


vector<complex<double>> idft(vector<complex<double>>&signal)
{
    vector<complex<double>>ans;
    int N = signal.size();
    cout << N << "\n\n";
    complex<double> i = complex<double>(0,1);
    // let zeta = e^(-2*pi*i/N);
    complex<double>zeta = exp(2.0*pi/N*i);

    for(int k=0; k<N; k++)
    {
        complex<double>sum = 0;
        for(int n=0; n<N; n++)
        {
            sum += multiplycomplex(signal[n], pow(zeta, n*k));
        }
        // rounding the real and the imaginary parts of the answer if they are so close to zero (compared using eps)
        sum = complex<double>((AsmallerthanEPS(real(sum))?0:real(sum))/N, (AsmallerthanEPS(imag(sum))?0:imag(sum))/N);
        
        ans.push_back(sum);
    }
    return ans;
}


int main()
{
    int n;      cin>>n;
    vector<complex<double>>signal(n);
    for(int i=0; i<n; i++)
    {
        double real,imag;
        cin>>real>>imag;
        signal[i] = complex<double>(real,imag);
    }
    auto ans = (idft(signal));
    cout(ans);
}