#include <bits/stdc++.h>
#include "AudioFile.h"

using namespace std;

#define cin(vec) for (auto &i : vec) cin >> i
#define cout(vec) for (auto &i : vec) cout << i << " ";cout << '\n';
#define pi (3.14159265358979323846)
#define EPS 0.01

// Construct Signals functions
vector<double> construct_vector(double initial_value, double final_value, double step);
vector<double> construct_sin_signal(double amplitude, double frequency, double phase_rad, vector<double>&time);
vector<double> construct_cos_signal(double amplitude, double frequency, double phase_rad, vector<double>&time);
vector<double> complex_to_amplitude_spectrum(vector<complex<double>>& signal);
vector<double> add_signals(vector<double>& signal1, vector<double>& signal2);

// Window functions
vector<double> triang(vector<double>& signal);
vector<double> hamming(vector<double>& signal);
vector<double> hanning(vector<double>& signal);

// Prepare files to plot using gnuplot
void write_to_file(vector<double>& signal, string file_name);
void write_to_file(vector<double>& x_axis, vector<double>& y_axis, string file_name);       // to print x axis values.


// DSP Operations / functions
bool AsmallerthanEPS (double a);
complex<double>multiplycomplex(complex<double>a, complex<double>b);

vector<double> convolution(vector<double>&source_signal, vector<double>&impulse_response);

vector<complex<double>> dft(vector<double>&signal);
vector<complex<double>> idft(vector<complex<double>>&signal);
void fftshift(vector<double>& signal);
/////////////////////////////////////////////////////////////////

// Frequency Response 
vector<double> complex_to_abs(vector<complex<double>>& signal);
vector<pair<complex<double>, double>>freqz(vector<double>B, vector<double>A, double N);
vector<double> to_db(vector<double>&h_amp);
/////////////////////////////////////////////////////////////////

//Filtering Functions
vector<double> to_db(vector<double>&h_amp);
vector<double> construct_filter_coeffiecients(vector<double>& pos_values, double zero_value);
vector<double> fir_filter(vector<double>& B, vector<double>& signal);
/////////////////////////////////////////////////////////////////

vector<complex<double>> fft(vector<complex<double>>& p);
vector<complex<double>> double_to_complex(vector<double>& p);

AudioFile<double> audioFile;

int main()
{
    audioFile.load ("C:/Users/Mohamed Amr/Downloads/Music/IFeelGood_10k.wav");
    audioFile.printSummary();
    int channel = 0;
    int numSamples = audioFile.getNumSamplesPerChannel();
    double fs = audioFile.getSampleRate();
    vector<double>signal(numSamples);
    for (int i = 0; i < numSamples; i++)
    {
        signal[i] = audioFile.samples[channel][i];
    }
    cout << "signal length = " << signal.size() << "\n";
    vector<double>t = construct_vector(0, 45, 1.0/fs);
    write_to_file(t,signal, "signal.dat");
    vector<complex<double>>signal_c = double_to_complex(signal);
    vector<complex<double>>dft_signal = fft(signal_c);
    int n = signal_c.size();
    vector<double> dft_signal_abs = complex_to_amplitude_spectrum(dft_signal);
    vector<double> frequency_axis = construct_vector(-fs/2, fs/2, fs/(n));
    fftshift(frequency_axis);
    write_to_file(frequency_axis, dft_signal_abs, "signal_dft.dat");
    double fc = 2000;                                    // Hz
    double omega_c = 2.0 * pi * fc / fs;
    int number_taps = 3301;
    int M = (number_taps-1)/2;
    double h_zero = omega_c / pi;
    vector<double>pos;                                  // pos values
    for(int i=1; i<=M; i++) 
        pos.push_back(sin(omega_c * i)/(i*pi));

    vector<double>B = construct_filter_coeffiecients(pos, h_zero);
    vector<double>A = {1};
    B = hamming(B);
    vector<pair<complex<double>, double>> output = freqz(B, A, 512);
    vector<complex<double>> H;                  // complex frequency Response 
    vector<double> W;                           // corresponding angle in radians
    for(auto& i:output)
    {
        H.push_back(i.first);
        W.push_back(i.second);
    }
    vector<double>H_amp = complex_to_abs(H);
    vector<double>h_amp_db = to_db(H_amp);
    for(int i=0; i<W.size(); i++)               // convert angle(rad) to frequency in Hz
    {
        W[i] *= fs / 2 / pi;
    }
    write_to_file(W, h_amp_db, "frequency response.dat");
    vector<double> filtered_signal = fir_filter(B, signal);
    write_to_file(t, filtered_signal, "Filtered Signal.dat");
    vector<complex<double>> filtered_signal_c = double_to_complex(filtered_signal);
    filtered_signal_c = fft(filtered_signal_c);
    filtered_signal = complex_to_amplitude_spectrum(filtered_signal_c);
    write_to_file(frequency_axis, filtered_signal, "Filtered Signal dft.dat");

    // //write_to_file(t, signal, "four.dat");
    // vector<complex<double>>dft_fsignal = dft(filtered_signal);
    // vector<double> dft_fsignal_abs = complex_to_amplitude_spectrum(dft_fsignal);
    // write_to_file(frequency_axis, dft_fsignal_abs, "three.dat");
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

vector<double> complex_to_amplitude_spectrum(vector<complex<double>>& signal)
{
    int N = signal.size();
    vector<double>ans(N);
    for(int i=0; i<signal.size(); i++)
    {
        ans[i] = abs(signal[i])*1/N;
    }
    return ans;
}

void write_to_file(vector<double>& x_axis, vector<double>& y_axis, string file_name)
{
    int N = min(x_axis.size(), y_axis.size());
    FILE* file_ptr = fopen(file_name.c_str(), "w");
    for(int i=0; i<N; i++)
    {
        fprintf(file_ptr, "\n%f \t %f", x_axis[i], y_axis[i]);
    }
    fclose(file_ptr);
}
void fftshift(vector<double>& signal)
{
    int n = signal.size();
    int mid = n/2;
    reverse(signal.begin(), signal.end());
    reverse(signal.begin(), signal.begin()+mid);
    reverse(signal.begin()+mid, signal.end());
}

vector<double> add_signals(vector<double>& signal1, vector<double>& signal2)
{
    int n = min(signal1.size(), signal2.size());
    vector<double>ans(n);
    for(int i=0; i<n; i++)
    {
        ans[i] = signal1[i] + signal2[i];
    }
    return ans;
}

vector<double> complex_to_abs(vector<complex<double>>& signal)
{
    int N = signal.size();
    vector<double>ans(N);
    for(int i=0; i<signal.size(); i++)
    {
        ans[i] = abs(signal[i]);
    }
    return ans;
}

vector<pair<complex<double>, double>>freqz(vector<double>B, vector<double>A, double N)
{
    vector<pair<complex<double>,double>>ans;
    vector<complex<double>>numerator(B.size());
    vector<complex<double>>denominator(A.size());
    complex<double> j = complex<double>(0, 1);
    complex<double>num;
    complex<double>denom;
    vector<double> list_of_omegas = construct_vector(0, 1*pi, pi/(N));
    for(auto& omega: list_of_omegas)
    {
        num = 0.0;
        denom = 0.0;
        for(double i=0; i<B.size(); i++)
        {
            num += B[i] * exp(-j*omega*i);
        }
        for(double i=0; i<A.size(); i++)
        {
            denom += A[i] * exp(-j*omega*i);
        }
        complex<double> temp = num / denom;
        ans.push_back(make_pair(temp, omega));
    }
    return ans;
}

vector<double> to_db(vector<double>&h_amp)
{
    vector<double>temp;
    for(auto& i:h_amp)
    {
        double h = 20 * log10(i);
        temp.push_back(h);
    }
    return temp;
}

vector<double> construct_filter_coeffiecients(vector<double>& pos_values, double zero_value)
{
    vector<double>first = pos_values;            // pos values
    vector<double>second = first;               // neg values
    reverse(second.begin(), second.end());
    second.push_back(zero_value);
    first.insert(first.begin(), second.begin(), second.end());
    return first;
}

vector<double> fir_filter(vector<double>& B, vector<double>& signal)
{
    int number_coefficients = B.size();
    vector<double>ans(signal.size(), 0.0);
    vector<double> buffer(number_coefficients, 0);          // initialize buffer with zeros

    for(int j=0; j<signal.size(); j++){
        double temp = 0.0;
        for(int i=number_coefficients-1; i>0; i--)
        {
            buffer[i] = buffer[i-1];
        }
        buffer[0] = signal[j];

        for(int i=0; i<number_coefficients; i++)
        {
            temp += B[i] * buffer[i];
        }
        ans[j] = temp;
    }
    return ans;
}

vector<complex<double>> fft(vector<complex<double>>& p)
{
    int n = p.size();
    if(n==1){
        return p;
    }
    int nn = pow(2, (int)(ceil(log2(p.size()))));   // largest power of two smaller than or equal to length of p 
    p.resize(nn);
    n = p.size();
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

vector<complex<double>> double_to_complex(vector<double>& p)
{
    vector<complex<double>>temp;
    for(auto& i:p)
        temp.push_back(complex<double>(i, 0));
    return temp;
}