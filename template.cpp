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
vector<double> construct_sin_signal(double amplitude, double frequency, double phase_rad, vector<double>&time);
vector<double> construct_cos_signal(double amplitude, double frequency, double phase_rad, vector<double>&time);
vector<double> convolution(vector<double>&source_signal, vector<double>&impulse_response);

int main()
{
    vector<double> temp = construct_vector(0, 1, 1/100.0);
    vector<double> temp2 = construct_sin_signal(1, 10, 0, temp);
    write_to_file(temp, "one.dat");
    write_to_file(temp2, "two.dat");
    vector<double> temp3 = hamming(temp2), temp4 = triang(temp2);
    write_to_file(temp3, "three.dat");
    write_to_file(temp4, "four.dat");
    cout(temp2);
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
vector<double> construct_sin_signal(double amplitude, double frequency, double phase_rad, vector<double>&time_vector)
{
    vector<double>temp(time_vector.size());
    double w = 2 * pi * frequency;
    for(int i=0; i<time_vector.size(); i++)
    {
        temp[i] = amplitude * sin(w*time_vector[i] + phase_rad);
    }
    return temp;
}
vector<double> construct_cos_signal(double amplitude, double frequency, double phase_rad, vector<double>&time_vector)
{
    vector<double>temp(time_vector.size());
    double w = 2 * pi * frequency;
    for(int i=0; i<time_vector.size(); i++)
    {
        temp[i] = amplitude * cos(w*time_vector[i] + phase_rad);
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

vector<double> convolution(vector<double>&source_signal, vector<double>&impulse_response)
{
    int n = source_signal.size();
    int m = impulse_response.size();
    vector<double> output_signal(n+m, 0);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            output_signal[i] += source_signal[i] * impulse_response[j];
        }
    }
    return output_signal;
}
