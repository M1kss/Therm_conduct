#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

/*static int N = 1000;
static int M = 1000;
static int K = 1000;
static float a = 0.5;
static float b = 0.5;*/

//static int T = 5; //

struct Data
{
    float a;
    float b;

    double tau;
    double hx;
    double hy;
    unsigned int N;
    unsigned int M;
    unsigned int K;
};

class DataMaintainer
{
public:

    double** f; //K*M
    double** mu_right; //M*N
    double** mu_top; //K*N
    double** mu_bottom; //K*N

    double** layer; //K*M
    double** layerSemi; //K*M

    Data* data;


    explicit DataMaintainer(const string & file)
    {
        data = new Data;

        ifstream in;
        in.open(file);
        in >> data->N;
        cout << data->N << " ";
        in >> data->K;
        cout << data->K << " ";
        in >> data->M;
        cout << data->M << " ";
        in >> data->a;
        in >> data->b;
        in >> data->tau;
        in >> data->hx;
        in >> data->hy;

        f = new double*[data->K+1];
        for (int k=0; k<=data->K; k++)
        {
            f[k] = new double[data->M+1];
        }

        mu_right = new double*[data->M+1];
        for (int m=0; m<=data->M; m++)
        {
            mu_right[m] = new double[data->N+1];
        }

        mu_top = new double*[data->K+1];
        for (int k=0; k<=data->K; k++)
        {
            mu_top[k] = new double[data->N+1];
        }

        mu_bottom = new double*[data->K+1];
        for (int k=0; k<=data->K; k++)
        {
            mu_bottom[k] = new double[data->N+1];
        }

        layer = new double*[data->K+1];
        for (int k=0; k<=data->K; k++)
        {
            layer[k] = new double[data->M+1];
        }

        layerSemi = new double*[data->K+1];
        for (int k=0; k<=data->K; k++)
        {
            layerSemi[k] = new double[data->M+1];
        }





        for (int k=0; k<=data->K; k++)
        {
            for(int m=0; m<=data->M; m++)
            {
                in >> layer[k][m];
            }
        }

        for (int m=0; m<=data->M; m++)
        {
            for(int n=0; n<=data->N; n++)
            {
                in >> mu_right[m][n];
            }
        }

        for (int k=0; k<=data->K; k++)
        {
            for(int n=0; n<=data->N; n++)
            {
                in >> mu_top[k][n];
            }
        }

        for (int k=0; k<=data->K; k++)
        {
            for(int n=0; n<=data->N; n++)
            {
                in >> mu_bottom[k][n];
            }
        }

        for (int k=0; k<=data->K; k++) //  4  ;  O?  @  >  A  B  >  B  K f=  57  0  2  8  A  8  B>  B n
        {
            for (int m = 0; m <= data->M; m++)
            {
                in >> f[k][m];
            }
        }

        in.close();
    }

    ~DataMaintainer()
    {
        for (int k=0; k<=data->K; k++)
        {
            delete [] f[k]; //K*M

            delete [] mu_top[k]; //K*N
            delete [] mu_bottom[k]; //K*N

            delete [] layer[k]; //K*M
            delete [] layerSemi[k]; //K*M
        }
        delete [] f; //K*M
        delete [] mu_top; //K*N
        delete [] mu_bottom; //K*N
        delete [] layer; //K*M
        delete [] layerSemi; //K*M

        for (int m=0; m<=data->M; m++)
        {
            delete[] mu_right[m]; //M*N
        }
        delete[] mu_right; //M*N

        delete data;
    }

};


class MathHandler
{
public:
    static void calculate_matrixSemi(int n, int m, int m_0, double** recvSemi, double** matrixSemi, Data* data, double** f, double** mu_right)
    {
        matrixSemi[0][0] = 0;
        matrixSemi[1][0] = 1;
        matrixSemi[2][0] = 0;
        matrixSemi[3][0] = recvSemi[0][m+1] + data->tau/2 * (4/(data->hx * data->hx)*(recvSemi[1][m+1] - recvSemi[0][m+1]) + f[0][m]);
        for (int k=1; k<data->K; k++)
        {
            matrixSemi[0][k] = -data->a * (1 - 1/(2*double(k))) / (data->hx * data->hx);
            matrixSemi[1][k] = 2*(1/data->tau + data->a / (data->hx * data->hx));
            matrixSemi[2][k] = -data->a * (1 + 1/(2*double(k))) / (data->hx * data->hx);
            matrixSemi[3][k] = 2/data->tau*recvSemi[k][m+1] + (1 - data->a)/(data->hx * data->hx)*( (1 - 1/(2*double(k)))*recvSemi[k-1][m+1] - 2*recvSemi[k][m+1] + (1 + 1/(2*double(k)))*recvSemi[k+1][m+1] ) + 1 / (data->hy * data->hy) * ( recvSemi[k][m] - 2*recvSemi[k][m+1] + recvSemi[k][m+2] ) + 0.5*f[k][m_0 + m];
        }
        matrixSemi[0][data->K] = 0;
        matrixSemi[1][data->K] = 1;
        matrixSemi[2][data->K] = 0;
        matrixSemi[3][data->K] = mu_right[m][n];
    }

    static void calculate_vectorSemi(int n, int m, int m_0, int m_1, double** recvSemi, double** matrixSemi, double* vectorSemi, Data* data, double** f, double** mu_top, double** mu_bottom, double** mu_right)
    {
        bool top = (m == data->M);
        bool bottom = (m == 0);
        if (top)
        {
            for (int k=0; k<=data->K; k++)
            {
                vectorSemi[k] = mu_top[k][n];
            }
        } else
        {
            if (bottom)
            {
                for (int k=0; k<=data->K; k++)
                {
                    vectorSemi[k] = mu_bottom[k][n];
                }
            } else
            {
                calculate_matrixSemi(n, m, m_0, recvSemi, matrixSemi, data, f, mu_right);
                shuttle(vectorSemi, matrixSemi, matrixSemi[3], data->K);
                /*for (int k=0; k<=data->K; k++)
                {
                    vectorSemi[k] = matrixSemi[3][k];
                }*/
            }
        }
    }

    static void calculate_stripSemi(int n, int m_0, int m_1, double** recvSemi, double** sendSemi, double** matrixSemi, Data* data, double** f, double** mu_top, double** mu_bottom, double** mu_right)
    {
        for (int m=0; m<=m_1-m_0; m++)
        {
            calculate_vectorSemi(n, m, m_0, m_1, recvSemi, matrixSemi, sendSemi[m], data, f, mu_top, mu_bottom, mu_right);
        }
    }

    static void calculate_matrix(int n, int k, int k_0, double** recv, double** matrix, Data* data, double** f, double** mu_top, double** mu_bottom)
    {
        matrix[0][0] = 0;
        matrix[1][0] = 1;
        matrix[2][0] = 0;
        matrix[3][0] = mu_bottom[k][n];
        for (int m=1; m<data->M; m++)
        {
            matrix[0][m] = -data->b / (data->hy * data->hy);
            matrix[1][m] = 2*(1/data->tau + data->b / (data->hy * data->hy));
            matrix[2][m] = -data->b / (data->hy * data->hy);
            matrix[3][m] = (2/data->tau)*recv[k+1][m] + 1 /(data->hx * data->hx)*( (1 - 1/(2*double(k)))*recv[k][m] - 2*recv[k+1][m] + (1 + 1/(2*double(k)))*recv[k+2][m] ) + (1 - data->b) / (data->hy * data->hy) * ( recv[k+1][m-1] - 2*recv[k+1][m] + recv[k+1][m+1] ) + 0.5*f[k_0 + k][m];
        }
        matrix[0][data->M] = 0;
        matrix[1][data->M] = 1;
        matrix[2][data->M] = 0;
        matrix[3][data->M] = mu_top[k][n];
    }

    static void calculate_vector(int n, int k, int k_0, int k_1, double** recv, double** matrix, double* vector, Data* data, double** f, double** mu_right, double** mu_top, double** mu_bottom)
    {
        bool right = (k == data->K);
        bool left = (k == 0);
        if (right)
        {
            for (int m=0; m<=data->M; m++)
            {
                vector[m] = mu_right[m][n];
            }
        } else
        {
            if (left)
            {
                for (int m=0; m<=data->M; m++)
                {
                    vector[m] = recv[1][m] + data->tau/2 * (4/(data->hx * data->hx)*(recv[2][m] - recv[1][m]) + f[0][m]);
                }
            } else
            {
                calculate_matrix(n, k, k_0, recv, matrix, data, f, mu_top, mu_bottom);
                shuttle(vector, matrix, matrix[3], data->M);
            }
        }
    }

    static void calculate_strip(int n, int k_0, int k_1, double** recv, double** send, double** matrix, Data* data, double** f, double** mu_top, double** mu_bottom, double** mu_right)
    {
        for (int k=0; k<=k_1-k_0; k++)
        {
            calculate_vector(n, k, k_0, k_1, recv, matrix, send[k], data, f, mu_right, mu_top, mu_bottom);
        }
    }

    static void shuttle(double* X, double**  A, double* F, int l)
    {
        for(int i = 1; i < l+1; i++){
            A[1][i] -= A[0][i]*A[2][i-1]/A[1][i-1];
            F[i] -= A[0][i]*F[i-1]/A[1][i-1];
        }
        X[l] = F[l]/A[1][l];
        for(int j = l-1; j > -1; j--){
            X[j] = (F[j] - A[2][j]*X[j+1])/A[1][j];
        }
    }
};


int main() {
    DataMaintainer DM = DataMaintainer("C:\\Users\\Sashok\\CLionProjects\\untitled\\init.txt");

    auto recvSemi = new double *[DM.data->K + 1];
    for (int k = 0; k <= DM.data->K; k++) {
        recvSemi[k] = new double[DM.data->M + 3];
    }

    auto recv = new double *[DM.data->K + 3];
    for (int k = 0; k <= DM.data->K + 2; k++) {
        recv[k] = new double[DM.data->M + 1];
    }

    auto sendSemi = new double *[DM.data->M + 1];
    for (int m = 0; m <= DM.data->M; m++) {
        sendSemi[m] = new double[DM.data->K + 1];
    }

    auto send = new double *[DM.data->K + 1];
    for (int k = 0; k <= DM.data->K; k++) {
        send[k] = new double[DM.data->M + 1];
    }

    auto matrixSemi = new double *[4];
    for (int i = 0; i <= 3; i++) {
        matrixSemi[i] = new double[DM.data->K + 1];
    }

    auto matrix = new double *[4];
    for (int i = 0; i <= 3; i++) {
        matrix[i] = new double[DM.data->M + 1];
    }

    for (int k = 0; k <= DM.data->K; k++) {
        for (int m = 0; m <= DM.data->M; m++) {
            recvSemi[k][m + 1] = DM.layer[k][m];
        }
    }






    for (int n=0; n<DM.data->N; n++)
    {
        cout << "n = " << n << endl;
        MathHandler::calculate_stripSemi(n, 0, DM.data->M, recvSemi, sendSemi, matrixSemi, DM.data, DM.f, DM.mu_top, DM.mu_bottom, DM.mu_right);
        for (int k=0; k<=DM.data->K; k++)
        {
            for (int m=0; m<=DM.data->M; m++)
            {
                recv[k+1][m] = sendSemi[m][k];
            }
        }

        /*cout << endl;
        for (int k=0; k<=DM.data->K; k++)
        {
            for (int m=0; m<=DM.data->M; m++)
            {
                cout << sendSemi[m][k] << " ";
            }
            cout << endl;
        }*/

        //cout << "n = 1/2 + " << n << endl;

        MathHandler::calculate_strip(n, 0, DM.data->K, recv, send, matrix, DM.data, DM.f, DM.mu_top, DM.mu_bottom, DM. mu_right);
        for (int k=0; k<=DM.data->K; k++)
        {
            for (int m=0; m<=DM.data->M; m++)
            {
                recvSemi[k][m+1] = send[k][m];
            }
        }

        /*cout << endl;
        for (int k=0; k<=DM.data->K; k++)
        {
            for (int m=0; m<=DM.data->M; m++)
            {
                cout << send[k][m] << " ";
            }
            cout << endl;
        }*/
    }

    ofstream out;
    out.open("C:\\Users\\Sashok\\CLionProjects\\untitled\\out.txt");
    cout << endl;

    out <<"{";
    for (int k=0; k<=DM.data->K; k++)
    {
        for (int m=0; m<=DM.data->M; m++)
        {
            if(m==DM.data->M)
            {
                cout << send[k][m] << endl;
                //out << send[k][m] << endl;
            }else{
                cout << send[k][m] << " ";
                //out << send[k][m] << " ";
            }
            out << "{" << DM.data->hx*k << ", " << DM.data->hy*m << ", " << send[k][m] << "}";
            //out << DM.data->hx*k << "," << DM.data->hy*m << "," << send[k][m] << endl;
            if(k!=DM.data->K || m!=DM.data->M) out << ", ";
        }

    }

    out << "}";



    for (int k=0; k<=DM.data->K+2; k++)
    {
        delete [] recv[k];
    }

    delete [] recv;


    for (int k=0; k<=DM.data->K; k++)
    {
        delete [] recvSemi[k];
        delete [] send[k];
    }

    for (int m=0; m<=DM.data->M; m++)
    {
        delete [] sendSemi[m];
    }

    delete [] recvSemi;
    delete [] send;
    delete [] sendSemi;

    for (int k=0; k<=3; k++)
    {
        delete [] matrix[k];
        delete [] matrixSemi[k];
    }

    delete [] matrix;
    delete [] matrixSemi;

    return 0;
}
