#include <bits/stdc++.h>
using namespace std;

#define cin(vec) for (auto &i : vec) cin >> i
#define cout(vec) for (auto &i : vec) cout << i << " ";cout << '\n';
#define pi (3.14159265358979323846)

vector<double> construct_vector(double initial_value, double final_value, double step);
vector<double> triang(vector<double>& signal);
vector<double> hamming(vector<double>& signal);
vector<double> hanning(vector<double>& signal);
void write_to_file(vector<double>& signal, string file_name);

int main()
{
    vector<double> temp = construct_vector(0, 1, 0.1);
    write_to_file(temp, "one.dat");
    temp = triang(temp);
    write_to_file(temp, "two.dat");
    cout(temp);
    return 0;
}

vector<double> construct_vector(double initial_value, double final_value, double step)
{
    vector<double>temp;
    for(double i=initial_value; i<=final_value; i+=step)
    {
        temp.push_back(i);
    }
    return temp;    
}

vector<double> triang(vector<double>& signal)
{
    int n = signal.size();
    vector<double>w(n);
    vector<double>temp(n);

    for(int i=0; i<n; i++)
    {
        w[i] = 1 - abs(2*i-n+1)/(double)(n-1);
    }

    for(int i=0; i<n; i++)
    {
        temp[i] = signal[i] * w[i];
    }
    return temp;
}
vector<double> hamming(vector<double>& signal)
{
    int n = signal.size();
    vector<double>w(n);
    vector<double>temp(n);

    for(int i=0; i<n; i++)
    {
        w[i] = 0.54 - 0.46 * cos(2 * pi * i / (n-1));
    }

    for(int i=0; i<n; i++)
    {
        temp[i] = signal[i] * w[i];
    }
    return temp;
}
vector<double> hanning(vector<double>& signal)
{
    int n = signal.size();
    vector<double>w(n);
    vector<double>temp(n);

    for(int i=0; i<n; i++)
    {
        w[i] = 0.5 - 0.5 * cos(2 * pi * i / (n-1));
    }

    for(int i=0; i<n; i++)
    {
        temp[i] = signal[i] * w[i];
    }
    return temp;
}

void write_to_file(vector<double>& signal, string file_name)
{
    FILE* file_ptr = fopen(file_name.c_str(), "w");
    for(auto& sample:signal)
    {
        fprintf(file_ptr, "\n%f", sample);
    }
    fclose(file_ptr);
}