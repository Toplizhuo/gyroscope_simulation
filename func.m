% Script written by:
% Zhuo Li (zhuol7@student.unimelb.edu.au)
% The University of Melbourne

function X_sol=func(t,X) 

a_0 = X(1);
b_0 = X(2);
c_0 = X(3);
d_0 = X(4);
a_1 = X(5);
b_1 = X(6);
c_1 = X(7);
d_1 = X(8);

expression = [(7078141215445656726073042565750*b_1*c_1 + 4245518245688925142712718458880*b_1*d_1 + 2500195204575568522761652178811*b_1*c_1*cos(c_0)^2 - 130176979546908692175566929920*b_1*d_1*cos(c_0)^2 + 147638083790227500439682678784*cos(c_0)*sin(b_0)*sin(c_0) - 81539949657416952006196871655354*a_1*b_1*cos(b_0) + 2500195204575568522761652178811*a_1*b_1*cos(b_0)*cos(c_0)^2 + 2500195204575568522761652178811*a_1^2*cos(b_0)*cos(c_0)*sin(b_0)*sin(c_0) + 2500195204575568522761652178811*a_1*c_1*cos(c_0)*sin(b_0)*sin(c_0) - 130176979546908692175566929920*a_1*d_1*cos(c_0)*sin(b_0)*sin(c_0))/(44309045436431304366134957110552*sin(b_0));
              (178043574288384*sin(b_0))/1690241007126979 + (10598877589516231*a_1^2*sin(2*b_0))/27043856114031664 + (18454760473778437554960334848*cos(c_0)^2*sin(b_0))/5538630679553913045766869638819 - (2500195204575568522761652178811*b_1*c_1*sin(2*c_0))/88618090872862608732269914221104 + (8136061221681793260972933120*b_1*d_1*sin(2*c_0))/5538630679553913045766869638819 - (2923050467499601*a_1*c_1*sin(b_0))/13521928057015832 - (156986423377920*a_1*d_1*sin(b_0))/1690241007126979 + (2500195204575568522761652178811*a_1*c_1*cos(c_0)^2*sin(b_0))/44309045436431304366134957110552 - (16272122443363586521945866240*a_1*d_1*cos(c_0)^2*sin(b_0))/5538630679553913045766869638819 + (2500195204575568522761652178811*a_1^2*cos(b_0)*cos(c_0)^2*sin(b_0))/44309045436431304366134957110552 - (2500195204575568522761652178811*a_1*b_1*cos(b_0)*cos(c_0)*sin(c_0))/44309045436431304366134957110552;
              (73870227446580795977937019074899808721446395032*a_1*b_1 + 28262380361590614194780550292109443069369113530*a_1*b_1*cos(b_0)^2 - 36742143539608995028253775267975567227786146752*a_1*b_1*cos(c_0)^2 - 8865703545381186754902693514814798424291134750*b_1*c_1*cos(b_0) - 5317710542514931507672782311114922040900976640*b_1*d_1*cos(b_0) - 3131611649818579604886207711957873014596575183*b_1*c_1*cos(b_0)*cos(c_0)^2 + 163052766816461016901622055467362547411189760*b_1*d_1*cos(b_0)*cos(c_0)^2 + 33610531889790415423367567556017694213189571569*a_1*b_1*cos(b_0)^2*cos(c_0)^2 - 18371071769804497514126887633987783613893073376*a_1^2*cos(c_0)*sin(b_0)*sin(c_0) + 18371071769804497514126887633987783613893073376*b_1^2*cos(c_0)*sin(b_0)*sin(c_0) - 184923618087194692237333625726727459126116352*cos(b_0)*cos(c_0)*sin(b_0)*sin(c_0) + 15239460119985917909240679922029910599296498193*a_1^2*cos(b_0)^2*cos(c_0)*sin(b_0)*sin(c_0) - 3131611649818579604886207711957873014596575183*a_1*c_1*cos(b_0)*cos(c_0)*sin(b_0)*sin(c_0) + 163052766816461016901622055467362547411189760*a_1*d_1*cos(b_0)*cos(c_0)*sin(b_0)*sin(c_0))/(55499155676776298463810131440912025107553321656*sin(b_0));
              (207306110849994*a_1^2*sin(2*c_0))/1252546858776253 - (207306110849994*b_1^2*sin(2*c_0))/1252546858776253 - (414612221699988*a_1*b_1*sin(b_0))/1252546858776253 + (829224443399976*a_1*b_1*cos(c_0)^2*sin(b_0))/1252546858776253 - (414612221699988*a_1^2*cos(b_0)^2*cos(c_0)*sin(c_0))/1252546858776253];


X_sol = zeros(8,1);
 
X_sol(1) = a_1;
X_sol(2) = b_1;
X_sol(3) = c_1;
X_sol(4) = d_1;
X_sol(5) = expression(1);
X_sol(6) = expression(2);
X_sol(7) = expression(3);
X_sol(8) = expression(4);

X_sol = [X_sol(1); X_sol(2); X_sol(3); X_sol(4); X_sol(5); X_sol(6); X_sol(7); X_sol(8)];

end