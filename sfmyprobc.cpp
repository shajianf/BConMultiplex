
#include<iostream>
#include<time.h>
#include<stdlib.h>
#include<fstream>
#include<math.h>
#include<cstdlib>
#include<vector>
#include<random>

#include<list>
#include<stack>
#include<queue>

#define N 10000				//number of nodes in each layers.
#define M 10				//main degree for SF networks.

using namespace std;

ofstream Inf;
ofstream Act;
ofstream Fee;

vector <vector<int>> A(N);			//adjecence matrix for above layer (opinion layer) (threshold model).
vector <vector<int>> B(N);			//adjecance matrix for another layer (disease layer) (SIS model).

//double p = .007;				// for ER networks (if you like to generate the ER network, first active this line).

double p[N] = {0};    				//SF networks parameters.             		
int degree [N]={0};				//degree of each nodes.
//int sort_degree[N]={0};
//int sort_degreenum[N]={0};

double gama_Net = 2.7;				//SF network parameter.
double alpha = 1/(gama_Net-1);			//SF network parameter.

double N0 =1000;                    		//initial active nodes (for getting start to dynamics).

int ensemble =50;			
int tmax = 1000;                 		// time steps.

double beta = 0.0;               		// probability of a node get infected (for SIS model).
double muo = .80;                		// probability of recovery (for sIS model).

double gama = .20;               		//Immunazation coefficent (an active node get infected by small probabitity (gama*beta=beta_A (beta_A <= beta))).
double kapa = 0.50;              		//self-awareness coefficent (probability that an infected node get active in above layer).

vector <int> sis(N);             		//state of each node at disease layer; 1= infected , 0= suseptible at time t (disease layer, SIS model). 
vector <int> sis_up(N);          		//state of each node at time t+1.

vector <int> thresh(N);          		//state of each node at opinion layer; 1= active , 0= inactive at time t (opinion layer, threshold model). 
vector <int> thresh_up(N);			//state of each node at time t+1.

vector <double> fee (N);        		//adoption threshold (consider the special thrreshold value for every node). 
//double fee =0.20;
vector <double> BC (N);         		//betweenness centrality values.


int main(){

    clock_t runtime = clock();

    Inf.open("Inf-rho.txt");			//output the density of infected nodes.
    Act.open("act-rho.txt");			//output the density of active nodes.
    //Fee.open("degree .txt");
    srand(time(NULL));

        std::default_random_engine generator;			//using random library.
       // std::uniform_real_distribution<double> distribution(0.0,1.0);
        //std::normal_distribution<double> distribution(.50,.05);

      //  for(int g =0; g <N; g++){

           //   fee [g] = distribution(generator) ;

        //Fee << fee [g] <<endl;
        //}


        double count =0;

        while (beta < 1.01){
            //cout<<"++++++++++++++++++"<< " " << beta <<endl;
            //p = k/N;
            double beta_A = beta*gama;

            int ens[ensemble] ={0};			//for epidemic output.
            int ens_A[ensemble] ={0};			//for opinion output.
     
            for(int e=0 ; e <ensemble ; e++){  
            
                for(int b =0; b < N; b++){              //reset all states and adjecance matrixes.

                    A[b].clear();
                    B[b].clear();
                    sis [b]=0;
                    sis_up [b]=0;
                    thresh [b]=0;
                    thresh_up [b]=0;
                   // p[b]=0;
                   // degree[b]=0;
                   // sort_degree[b]=0;
                   // sort_degreenum[b]=0;
                   // sort_degree[2][b]=0;
                   // fee[b] = 0;
                  //  p[b]=0;
                }

                A.resize(N);
                B.resize(N);

              // ----------------------- ER networks generate---------------------------;
              /*  for (int i=0; i< N;i++){
                    for (int j=0; j<i; j++){

                        float ran = rand()/(1.+RAND_MAX);
                        if( ran < p){

                            A[i].push_back(j);
                            A[j].push_back(i);  

                        //degree[i]++;
                        //degree[j]++;
                        }

                        float ran2 = rand()/(1.+RAND_MAX);
                        if( ran2 < p ){

                            B[i].push_back(j);
                            B[j].push_back(i);
                        }
                    }
                }*/

		//------------------------- SF generate networks -------------------------------.
                double constant =0;                    
                double tot_prob = 0 ;

	            for (int i = 0; i < N ; i++){		//asign weight w(i) for node i according 
		            tot_prob += pow((i+1),-alpha);
	            }

	            // cout << "TOTAL_PROB = " << tot_prob << endl;

	            for (int i = 0; i < N ; i++){
		            p[i] = pow((i+1),-alpha) / double(tot_prob);
		        // constant += p[i];
	            }
	            // cout << "TOT = " << constant << endl;
                int u =0;

                while(u <2){

                    double mean_degree=0;
                    double total_degree = 0;

                    while (mean_degree != M ){

                        int x1[2] ={0};

			            for (int j=0; j<2 ; j++){

				            double rand_num1 = rand() / double(RAND_MAX);
				            long double sum_of_prob1 = 0.0000000;

				            for (int i = 0 ; i < N ; i++){	
				
					            sum_of_prob1 += p[i];
					            if (rand_num1 <= sum_of_prob1){
						            x1[j] = i ;
						            break;
					            }
				            }
			            }

                        bool flag0 =true;
                        if ( x1[0]==x1[1] )	flag0 =false;

			            bool flag1 = true; 

                        int s = x1[0];
                        int d = x1[1];

			            for (int i = 0; i < B[s].size(); i++){

				            if ( B[s][i]==x1[1] ){
					
                                flag1 = false;
                                break;
				            }
                
			            }

			            if (flag1 && flag0) {

                            if(u==0){
                        
			                    B[s].push_back(d);
			                    B[d].push_back(s);

                              //  degree[d]++;
                               // degree[s]++;

                            }else {

                                A[s].push_back(d);
                                A[d].push_back(s);

                              //  degree[d]++;
                               // degree[s]++;
                            }

			                total_degree +=2;
			                mean_degree = total_degree*1.0 / N * 1.0;

                        }
	                }//while 

                    u++;
                }// while u.

                // ---------------------------------- phi allocating (threshold values) -------------------------------.
                double fee_sort[N] = {0};

                std::uniform_real_distribution<double> distribution(0,1);
                //std::normal_distribution<double> distribution(.5,.05);
                for(int g=0; g<N; g++){

                    fee [g] = distribution(generator);

                  //  Fee<< g << "\t" << fee[g] <<endl;
                }

               // Fee << "++++++++++++++++++++"<<endl;

                 
              //-------------------------- sort phi values of each node from big to small ----------------------------.
		    
               for(int c =0; c<N ; c++){                  

                    fee_sort[c] = fee[c]; 
                    //sort_degreenum[c] = c;

                  //  Fee << c << "\t" <<sort_degree[c]<< endl;

                }

                for(int v=0; v< N ; v++){
                    for(int n=0; n<N-1 ;n++){
                        double x=0;
                       // int y =0;

                        if(fee_sort [n] < fee_sort [n+1]){

                            x = fee_sort [n];
                            //y = sort_degreenum[n];

                            fee_sort [n] = fee_sort [n+1];
                           // sort_degreenum [n] = sort_degreenum [n+1];

                            fee_sort [n+1] = x;
                           // sort_degreenum [n+1] = y;
                        }

                    }
                }

                //++++++++++++++++++++++++++++++++++++++++++++ compute betweenness centrality.+++++++++++++++++++++++++++++++++++++++++++++++++++++.
		//using Brandes algorithem (2001).
		    
                for (int u1 = 0; u1 < N ; u1++){
                    BC[u1] = 0;
                    // sigma[u1]=0;
                }

                for (int s=0; s<N;s++){

                    stack <int> myStack;          	//empty stack.
                    list <int>  pred[N];                //empty list.

                    vector <int> dist(N);
                    vector <double> sigma(N);

                    for (int t1=0; t1<N; t1++){

                    sigma [t1] = 0;
                    dist [t1] = -1;

                    }

                    sigma [s] = 1;
                    dist [s] = 0;

                    queue <int> myqueue ;
                    myqueue.push(s);

        

                    while (! myqueue.empty()){

                        int v =0;
                        v = myqueue.front();
                        myqueue.pop();

                        myStack.push(v);

                        for(int u2=0;u2 < A[v].size(); u2++){ //omega found for first time?

                            int omega =0;

                            omega = A[v][u2];

                            if (dist [omega] < 0){

                            myqueue.push(omega);

                            dist[omega] = dist[v]+1;

                            }
                            if (dist[omega] == 1 + dist[v] ){

                                sigma[omega] = sigma[omega]+ sigma[v];
                                pred[omega].push_back(v);
                     
                            }

                        }//end adj.size()    

                    }//while queue. 
        

                    vector <double> delta(N);

                    for (int u3=0; u3<N; u3++){
            
                        delta[u3]=0;
                    }
        
                    while (! myStack.empty()){

                    // cout << "slam"<<endl;
                        int omega1 = 0;
                        omega1 = myStack.top();
                        myStack.pop();

            
            
                        for(list <int>::iterator u4= pred[omega1].begin(); u4 != pred[omega1].end() ; u4++){

                            int v1 =0;
                            v1 = *u4;
                
                            delta[v1] = (delta[v1] + (sigma[v1]/sigma[omega1])*(1.0 + delta[omega1]));

                        }
            

                        if (omega1 != s){

                            BC[omega1] = BC[omega1] + delta [omega1];
                        }


                    }
        
                  //  cout << s<<endl;

                }//end for s.

                //for (int u5 =0; u5<N; u5++){

                   // betce << u5 << "\t" << BC[u5] <<endl;
                //}

                //++++++++++++++++++++++++++++++ end of BC ++++++++++++++++++++++ . 
		//---------------------------- sort nodes by their BC values ------------------------.
               int sort_degreenum[N]={0};
                double sort_degree[N]={0};
                    
                for(int c =0; c<N ; c++){                   

                    sort_degree[c] = BC[c]; 
                    sort_degreenum[c] = c;

                  //  Fee << c << "\t" <<sort_degree[c]<< endl;

                }

               // Fee << "++++++++++++++++++++"<<endl;

                for(int v=0; v< N ; v++){
                    for(int n=0; n<N-1 ;n++){
                        int x=0;
                        int y =0;

                        if(sort_degree[n] < sort_degree[n+1]){

                            x = sort_degree[n];
                            y = sort_degreenum[n];

                            sort_degree [n] = sort_degree [n+1];
                            sort_degreenum [n] = sort_degreenum [n+1];

                            sort_degree [n+1] = x;
                            sort_degreenum [n+1] = y;
                        }

                    }
                }

                /*for(int uu=0; uu<N;uu++){

                    Fee<<sort_degreenum[uu] <<"\t" <<sort_degree[uu] <<endl;
                }
            
                Fee << "#########*************######### "<< endl; */
               /* int y=0;
                while(y<N0 ){           //initilaze in opinion layer

                    int O=0;
                    O = rand()%N;
                   // O = sort_degreenum[y];
                    if (thresh_up[O]==0){

                        //sis[O]=1;
                        //sis_up[O]=1;
                        thresh[O]=1;
                        thresh_up[O]=1;
                       // Fee << O <<endl;
                       // Fee << y << "\t" << sort_degree[y] << "\t" << sort_degreenum[y]<<"\t" << degree[O] << "\t" << O<<endl;
                        y++;
                        
                    }
                   
                  //  Fee << y << "\t" << sort_degree[y] << "\t" << sort_degreenum[y]<<"\t" << degree[O] << "\t" << O<<endl;
                  //  y++;
                } */


                int y1=0;
                while(y1 <N0) {                 	//initilaze disease layer.

                    int o1=0;
                    o1 = rand()%N;
                  //  o1 = sort_degreenum[y1];

                    if (sis_up[o1]==0){

                        sis_up[o1]=1;
                        sis[o1] =1;

                        y1++;

                    }
                }
		
		//---------------------------------- sort phi values with node's betweenness centrality --------------.    
		    
                for(int qq =0; qq<N; qq++){              

                    double o2 =0;
                    o2 = sort_degreenum[qq];

                    fee[o2] = fee_sort[qq];

                  //  Fee << o2 << "\t" << fee[o2]<< "\t" << sort_degree[qq] <<endl;
                }



               // Fee << "ooooooOOOooooooo"<<endl;
		//--------------------------------- start dynamics --------------------------------.
                int t=0;  				
                while(t<tmax){
			//------------------------- threshold model ---------------------------------.
                    for (int f=0; f<N; f++){            

                       
                        int act_nei =0;
                        int k_num =0;

                    
                        k_num = A[f].size();

                        for (int f1 =0; f1 < A[f].size() ; f1++){           //calculate the active neighbours.

                            int q1=0;
                            q1 = A[f][f1];

                            if(thresh[q1]==1){

                                act_nei++;
                            }

                        }
                        //}
                        //Act << fee <<"\t" << ((double) act_nei) <<"\t " << k_num <<endl;
                        if(fee_sort[f]  < ((double) act_nei/k_num)){

                            thresh_up[f] =1;
                        }else{

                            thresh_up[f] =0;
                        }
			// ------------------- end threshold model -------------------------------.

                        if((sis[f]==1)&&(thresh[f]==0)){                //self-awerness with kapa prob.

                            float ran3 = rand()/(1.+RAND_MAX);

                            if (ran3 < kapa){

                                thresh_up [f]=1;

                               // cout << "oooohhhhh!!!!" <<endl; 
                            }
                        }

			//-------------------------------------- SIS model --------------------------.
                        float ran1= rand()/(1.+RAND_MAX);

                        if ((sis[f]==1) && (ran1 < muo)){           //get recovery with mu prob.
                            //cout <<"hello" << endl;
                            sis_up[f]=0;
                            //thresh_up[f]=0;
                        }  

                        for(int h=0;h< B[f].size(); h++){

                        

                            bool flag =false;
                            int q=0;
                            q = B[f][h];

                        

                            float ran2 = rand()/(1.+RAND_MAX);

                           if((thresh[f]==1) && (sis[f]==0) && ran2 < (beta_A * sis[q])){  //active node get infected with beta_A prob.
                                
                                sis_up[f]=1;
                                flag = true;

                            }  
                            if(flag) break; 

                        

                            if((thresh[f]==0) && (sis[f]==0) && (ran2 < (beta * sis[q]))){  //inactive node get infected with beta prob.

                                sis_up [f]=1;
                                flag = true;

                            }
                            if(flag) break;
                        }
                    
                    }
			//------------------------ end SIS model ---------------------------------.
                
                    for( int u=0; u<N ; u++){		//change steps for new step.
                    
                        sis[u] = sis_up[u];
                        thresh [u] = thresh_up[u];
                    }
                    t++;

                   /* if(t % (400) ==0){

                        for (int h3 =0; h3<N; h3++){

                            if((fee[h3]< 0.60)){

                                fee[h3] += 0.12;
                            }
                        }

                    }*/

                }//end of while t.
        
                int densityofinf =0;
                int densityofact =0;
                for(int b=0; b < N; b++){

                    densityofinf += sis_up[b];
                    densityofact += thresh_up[b];
                }

                ens[e] = densityofinf;
                ens_A[e] = densityofact;
            
            }//end of ensemble.

            double inf =0;
            double act =0;
            for(int a=0; a<ensemble ;a++){
                inf +=ens[a];
                act +=ens_A[a];

                ens[a]=0;
                ens_A[a]=0;
            }
            
            Inf /*<<kapa <<'\t'*/ << beta <<'\t'<< inf/double(N * ensemble)<< endl;		//txt output.
            Act /*<<kapa << '\t' */<< beta <<"\t" << act/double(N * ensemble)<< endl;		//txt output.
            //cout<< fee <<'\t' << beta << '\t' << act/double(N * ensemble ) <<"\t" << inf/double(N * ensemble ) <<  endl;
            cout << beta <<endl;

            beta += 0.01;
            //domineEnd += 0.02;
           // dominefirst -= 0.01;
        }

      //  Inf << endl;
       // Act << endl;

       // cout << kapa << endl;

       // kapa += .02;

    //}

    cout << "run time = " << (double) (clock() - runtime) / (CLOCKS_PER_SEC * 60) << " min" << endl;

 return 0;

}
