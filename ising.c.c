#include <math.h>
#include <stdio.h>
#include "gsl_rng.h" //Libreria para generación de números aleatorios
#define N 64
#define pasos 1000001
#define fps 100

gsl_rng *tau;

int main()
{   
    int s[N][N];
    int i,j,k,l,n,m,o,sum,b;
    int semilla=76854;       //Semilla de generación de números
    double dE,P,x;
    extern gsl_rng *tau;    //Puntero al estado del número aleatorio
    double unpmontecarlo, pasosmontecarlo;
    int resto,contador;
    double T;
    double mag,summag,Emed,sumEmed,Emedcuadrado,sumEmedcuadrado,cespec,fcorrelacion,fcorrelacioni[N],sumfcorrelacion[N-1];
    unpmontecarlo=N*N;          //Definimos un paso montecarlo, que es el número de spines
    pasosmontecarlo=pasos;  //Definimos el número de pasos montecarlo
    FILE *f1,*f2,*f3,*f4,*f5;
    
    tau=gsl_rng_alloc(gsl_rng_taus);    //Inicializamos el puntero
    gsl_rng_set(tau,semilla);           //Inicializamos la semilla
    T=1.5;
    f1=fopen("datos.txt", "w");
    f2=fopen("mag.txt","w");
    f3=fopen("Emed.txt","w");
    f4=fopen("Emedcuadrado.txt","w");
    f5=fopen("fcorrelacion.txt","w");
    fprintf(f1, "Temp\t\tMagnet\t\tEmedia\t\tCespecífico\n");
for(o=0;o<10; o++)
{
    contador=0;
    summag=0.0;
    sumEmed=0.0;
    sumEmedcuadrado=0.0;
    fprintf(f2, "%lf\n\n", T);
    fprintf(f3, "%lf\n\n", T);
    fprintf(f4, "%lf\n\n", T);
    for(b=0;b<N-1; b++)
    {
        sumfcorrelacion[b]=0.0;
    }
    for(b=0;b<N-1; b++)
    {   
        fcorrelacioni[b]=0.0;
    }
    //Inicializamos los spines
    for (i=0;i<N; i++)
    {
        for (j=0;j<N; j++)
        {
            s[i][j]=1;
        }
    }
    //Comenzamos a iterar.
    //Queremos guardar los datos para cada paso montecarlo
    for(i=1;i<pasosmontecarlo; i++)
    {        
        for(j=0;j<unpmontecarlo; j++)
        {
            n=gsl_rng_uniform_int(tau,N);   //Asigno a n un valor entero entre [0,N-1]
            m=gsl_rng_uniform_int(tau,N);   //Igual para m
        
            //Calculamos la diferencia de energía
            //Imponemos las condiciones de contorno periódicas
            if(n!=0 && m!=0 && n!=N-1 && m!=N-1)
            {
                dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
            }
            //Con el if anterior quitamos todas las posiciones que no son bordes
            //Si estoy en el borde o en una de las esquinas entra en el siguiente else
            else
            {
                if(n==0 && m==0 || n==N-1 && m==N-1 || n==0 && m==N-1 || n==N-1 && m==0)
                //Con este if selecciono si estoy en una de las esquinas
                {
                    if(n==0 && m==0)
                    {
                        dE=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][N-1]);
                    }
                    if(n==N-1 && m==N-1)
                    {
                        dE=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
                    }
                    if(n==0 && m==N-1)
                    {
                        dE=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][0]+s[n][m-1]);
                    }
                    if(n==N-1 && m==0)
                    {
                        dE=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
                    }
                }
                //Si no estoy en ninguna esquina y ha entrado en este else, significa que
                //Estoy en un borde, sin ser esquina, por tanto solo una de las dos, n o m
                //Tiene que tenerse en cuenta
                else
                {
                    if(n==0)
                    {
                        dE=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][m-1]);
                    }
                    if(m==0)
                    {
                        dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
                    }
                    if(n==N-1)
                    {
                        dE=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
                    }
                    if(m==N-1)
                    {
                        dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
                    }
                }
            }
            if(dE<=0)
            {
                P=1;
            }
            else
            {
                P=exp(-(1.0*dE)/T);
            }
            x=gsl_rng_uniform(tau);  //número aleatorio real [0,1]
            if(x<P)
            {
                s[n][m]=-s[n][m];
            }
            else
            {
                s[n][m]=s[n][m];
            }
        }
        resto=i%fps;    //Calculamos los datos cada 100 pasos montecarlo
        if(resto==0)
        {
            //Calculamos primero la magnetización
            sum=0.0;
            for(k=0;k<N; k++)
            {
                for(l=0;l<N; l++)
                {
                    sum=sum+s[k][l];
                }
            }
            sum=sqrt(sum*sum);
            mag=(sum*1.0)/(N*N);
            summag=summag+mag;
            //Calculamos la energía media
            Emed=0.0;
            //Con el siguiente for contamos con todas las casillas que no están en el borde
            for(k=1;k<N-1; k++)
            {
                for(l=1;l<N-1; l++)
                {
                    Emed=Emed+s[k][l]*(s[k][l-1]+s[k][l+1]+s[k-1][l]+s[k+1][l]);
                }
            }
            //Ahora contamos los bordes
            for(k=1;k<N-1; k++)
            {
                Emed=Emed+s[0][k]*(s[0][k-1]+s[0][k+1]+s[N-1][k]+s[1][k]);  //Borde de arriba
                Emed=Emed+s[N-1][k]*(s[N-1][k-1]+s[N-1][k+1]+s[N-2][k]+s[0][k]);    //Borde de abajo
                Emed=Emed+s[k][0]*(s[k][N-1]+s[k][1]+s[k-1][0]+s[k+1][0]);  //Borde de la izquierda
                Emed=Emed+s[k][N-1]*(s[k][N-2]+s[k][0]+s[k-1][N-1]+s[k+1][N-1]);    //Borde de la derecha
            }
            //Finalmente contamos las esquinas            
            Emed=Emed+s[0][0]*(s[0][N-1]+s[0][1]+s[N-1][0]+s[1][0]);    //Superior izquierda
            Emed=Emed+s[0][N-1]+(s[0][N-2]+s[0][0]+s[N-1][N-1]+s[1][N-1]);  //Superior derecha
            Emed=Emed+s[N-1][0]*(s[N-1][N-1]+s[N-1][1]+s[N-2][0]+s[0][0]);  //Inferior izquierda
            Emed=Emed+s[N-1][N-1]*(s[N-1][N-2]+s[N-1][0]+s[N-2][N-1]+s[0][N-1]);    //Inferior derecha
            Emed=-(1.0/2.0)*Emed;
            sumEmed=sumEmed+Emed;
            Emedcuadrado=Emed*Emed;
            sumEmedcuadrado=sumEmedcuadrado+Emedcuadrado;;
            //Calculamos la función de correlacion
            for(b=1;b<N; b++)
            {
                fcorrelacion=0;
                for(k=0;k<N-b; k++)
                {
                    for(l=0;l<N; l++)
                    {
                        fcorrelacion=fcorrelacion+s[k][l]*s[k+b][l];
                    }
                }
                for(k=N-b;k<N; k++)
                {
                    for(l=0;l<N; l++)
                    {
                        fcorrelacion=fcorrelacion+s[k][l]*s[k-N+b][l];
                    }
                }
                sumfcorrelacion[b-1]=sumfcorrelacion[b-1]+fcorrelacion;
            }

            //Metemos los datos en el fichero para poder procesarlos luego
            fprintf(f2, "%lf\n", mag);
            fprintf(f3, "%lf\n", Emed);
            fprintf(f4, "%lf\n", Emedcuadrado);
            contador=contador+1;
        }
    }
    mag=summag/contador;
    Emed=(1.0*sumEmed)/(contador);
    Emedcuadrado=(1.0*sumEmedcuadrado)/contador;
    cespec=((Emedcuadrado-Emed*Emed)*(Emedcuadrado-Emed*Emed))/(N*N*T);
    Emed=(Emed*1.0)/(2*N);
    for(b=0;b<N-1; b++)
    {
        sumfcorrelacion[b]=sumfcorrelacion[b]/(contador*N*N);
    }
    fprintf(f2, "\n");
    fprintf(f3, "\n");
    fprintf(f4, "\n");
    fprintf(f1, "%lf\t%lf\t%lf\t%lf\n", T,mag,Emed,cespec);
    fprintf(f5, "%lf\n", T);
    for(b=1;b<N; b++)
    {
        fprintf(f5, "%i\t%lf\n", b,sumfcorrelacion[b-1]);
    }

    T=T+0.22;
}
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    return (0);
}
