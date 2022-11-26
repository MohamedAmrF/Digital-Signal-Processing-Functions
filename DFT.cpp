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

    // let zeta = e^(-2*pi*i/N);
    complex<double>zeta = exp(-2*pi/N*complex<double>(0,1));

    for(int i=0; i<N; i++)
    {
        complex<double>sum = 0;
        for(int j=0; j<N; j++)
        {
            complex<double>temp;
            temp = signal[j] * pow(zeta, j*i);
            sum += temp;
        }
        // rounding the real and the imaginary parts of the answer if they are so close to zero (compared using eps)
        sum = complex<double>((AsmallerthanEPS(real(sum))?0:real(sum)), (AsmallerthanEPS(imag(sum))?0:imag(sum)));
        
        ans.push_back(sum);
    }
    return ans;
}

void solve()
{       
    int n;      cin>>n;
    vector<double>signal(n);
    cin(signal);
    auto ans = (dft(signal));
    cout(ans);
}

int main()
{
    int n;      cin>>n;
    vector<double>signal(n);
    cin(signal);
    auto ans = (dft(signal));
    cout(ans);
}