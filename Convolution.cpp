#include<bits/stdc++.h>

using namespace std;
void write_to_file(vector<double>& signal, string file_name);
vector<double> convolution(vector<double>&source_signal, vector<double>&impulse_response);

int main()
{
    vector<double> signal = {0, 1, 2, 3};
    vector<double> impulse = {1, 2, 1};
    vector<double> result = convolution(signal, impulse);
    write_to_file(signal, "input.dat"); // to plot using gnu plot
    write_to_file(Impulse_response, "impulse.dat"); // to plot using gnu plot
    write_to_file(result, "output_signal.dat"); // to plot using gnu plot
    return 0;
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

void write_to_file(vector<double>& signal, string file_name)
{
    FILE* file_ptr = fopen(file_name.c_str(), "w");
    for(auto& sample:signal)
    {
        fprintf(file_ptr, "\n%f", sample);
    }
    fclose(file_ptr);
}