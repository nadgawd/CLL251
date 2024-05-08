#include<bits/stdc++.h>
using namespace std;

float interpolate(float x, float x0, float x1, float y0, float y1)
{
    return y0 + (x - x0)*(y1 - y0)/(x1 - x0);
}

int main()
{
    float W1, W2, L1, L2, b1, b2, Vfr1, Vfr2, Tint, Text, Tc, Th, v1, v2, K1, K2, Pr1, Pr2, No1, No2, At1, At2;
    cout<<"Enter Width(in m) of Internal and External Heat Sink: "<<endl;
    cin>>W1>>W2;
    cout<<"Enter Length(in m) of Internal and External Heat Sink: "<<endl;
    cin>>L1>>L2;
    cout<<"Enter Breadth(in m) of Internal and External Heat Sink: "<<endl;
    cin>>b1>>b2;
    cout<<"Enter Volumetric Flow Rates(in cms) of Internal and External Heat Sink: "<<endl;
    cin>>Vfr1>>Vfr2;
    cout<<"Enter Internal and External Air Temperatures: "<<endl;
    cin>>Tint>>Text;
    cout<<"Enter Cold and Hot Juntion Temperatures: "<<endl;
    cin>>Tc>>Th;

    // Cross Flow Rates
    float Ac1 = b1*W1;
    float Ac2 = b2*W2;
    // cout<<"Ac1 = "<<Ac1<<endl;
    // cout<<"Ac2 = "<<Ac2<<endl;

    // Velocities
    float U1 = Vfr1/(2*Ac1);
    float U2 = Vfr2/(2*Ac2);
    // cout<<"U1 = "<<U1<<endl;
    // cout<<"U2 = "<<U2<<endl;

    float theta1 = Tint - Tc;
    float theta2 = Th - Text;

    float data[35][4];
    // T v K Pr
    ifstream rfile("thermophysicalprops.txt");
    
    for   (int i = 0; i < 35; i++) {
        for (int j = 0; j < 4; j++) {
            rfile >> data[i][j];
        }
    }
    
    // for (int i = 0; i < 35; i++) {
    //     for (int j = 0; j < 4; j++) {
    //         cout << data[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    float Tavg1 = (Tint + Tc)/2 + 273.15;
    float Tavg2 = (Text + Th)/2 + 273.15;

    for (int i=0; i<35; i++){
        if(data[i][0]<Tavg1){
            continue;
        }
        else if(data[i][0]>Tavg1){
            v1 = interpolate(Tavg1, data[i-1][0], data[i+1][0], data[i-1][1], data[i][1]);
            K1 = interpolate(Tavg1, data[i-1][0], data[i+1][0], data[i-1][2], data[i][2]);
            Pr1 = interpolate(Tavg1, data[i-1][0], data[i+1][0], data[i-1][3], data[i][3]);
            break;
        }
    }
    for (int i=0; i<35; i++){
        if(data[i][0]<Tavg2){
            continue;
        }
        else if(data[i][0]>Tavg2){
            v2 = interpolate(Tavg2, data[i-1][0], data[i+1][0], data[i-1][1], data[i][1]);
            K2 = interpolate(Tavg2, data[i-1][0], data[i+1][0], data[i-1][2], data[i][2]);
            Pr2 = interpolate(Tavg2, data[i-1][0], data[i+1][0], data[i-1][3], data[i][3]);
            break;
        }
    }

    // cout<<v1<<"  "<<K1<<"    "<<Pr1<<endl;
    // cout<<v2<<"    "<<K2<<"    "<<Pr2<<endl;

    // Properties of aluminium
    float Rhoal = 2702;
    float kal = 177;

    float Lc1 = L1/2;
    float Lc2 = L2/2;
    

    // Reynolds number
    float Re1 = U1*Lc1/v1;
    float Re2 = U2*Lc2/v2;
    // cout<<"Re1 = "<<Re1<<endl;
    // cout<<"Re2 = "<<Re2<<endl;

    // optimal fin spacing
    float Zopt1 = Lc1*3.24*pow(Re1, -0.5)*pow(Pr1, -0.25);
    float Zopt2 = Lc2*3.24*pow(Re2, -0.5)*pow(Pr2, -0.25);
    // cout<<"Zopt1 = "<<Zopt1<<endl;
    // cout<<"Zopt2 = "<<Zopt2<<endl;

    // average heat transfer coefficients
    float h1 = K1*0.664*pow(Re1, 0.5)*pow(Pr1, 1/3)/Lc1;   
    float h2 = K2*0.664*pow(Re2, 0.5)*pow(Pr2, 1/3)/Lc2;
    // cout<<"h1 = "<<h1<<endl;
    // cout<<"h2 = "<<h2<<endl;

    float q1max = 0;
    float q2max = 0;

    float t1max = 0;
    float t2max = 0;

    float Rto1 = 0;
    float Rto2 = 0;

    ofstream fk;
    fk.open("output1.txt");

    for (float t1 = 1; t1 < 201; t1 ++){
        // cout<<"loop 1"<<endl;
        float nf1 = W1/(Zopt1+0.00001*t1);            // number of fins
        float B1 = b1*pow(2*h1/(kal*0.00001*t1), 0.5);
        float Nf1 = tanh(B1)/B1;             // single finn efficiency
        float Af1 = 2*(L1+0.00001*t1)*b1;           //single fin area
        At1 = nf1*(Af1+L1*Zopt1);        // total surface area
        No1 = 1 - nf1*(Af1/At1)*(1-Nf1);  //overall efficiency
        float qtot1 = No1*At1*h1*theta1;        //total heat transfers
        // cout<<qtot1;
        if(qtot1>q1max){
            q1max = qtot1;
            t1max = 0.00001*t1;
            Rto1 = 1/(No1*h1*At1);
        }
        fk<<qtot1<<setw(20)<<t1*0.01<<endl;
    }
    fk.close();

    ofstream tt;
    tt.open("output2.txt");
    for (float t2 = 1; t2 < 201; t2++){
        // cout<<"loop 2"<<endl;
        float nf2 = W2/(Zopt2+0.00001*t2);            // number of fins
        float B2 = b2*pow(2*h2/(kal*0.00001*t2), 0.5);
        float Nf2 = tanh(B2)/B2;             // single finn efficiency
        float Af2 = 2*(L2+0.00001*t2)*b2;           //single fin area
        At2 = nf2*(Af2+L2*Zopt2);        // total surface area
        No2 = 1 - nf2*(Af2/At2)*(1-Nf2);  //overall efficiency
        float qtot2 = No2*At2*h2*theta2;        //total heat transfers
        if(qtot2>q2max){
            q2max = qtot2;
            t2max = 0.00001*t2;
            Rto2 = 1/(No2*h2*At2);
        }
        tt<<qtot2<<setw(20)<<t2*0.01<<endl;
    }
    tt.close();

    cout<<"q1 max = "<<q1max<<endl;
    cout<<"q2 max = "<<q2max<<endl;
    cout<<"t1 max = "<<t1max*1000<<" mm"<<endl;
    cout<<"t2 max = "<<t2max*1000<<" mm"<<endl;
    cout<<"R1 max = "<<Rto1<<" W"<<endl;
    cout<<"R2 max = "<<Rto2<<" W"<<endl;
    rfile.close();
}