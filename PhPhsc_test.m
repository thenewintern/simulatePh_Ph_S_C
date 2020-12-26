function []=PhPhsc_test(t,prob)

global m_a m_b total_number_of_equations number_of_equations_per_phase c s;


dummy_prob=prob(:,1:total_number_of_equations);

%@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 1 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%%%%%%%%%%%%%%% Test for E'[N_i^0(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0_subspace1_kfe=0;

for l=1:m_a
    index=(l-1)*number_of_equations_per_phase;
    total_done=0;
    
    for n=0:s-1
        number_of_equations=nchoosek(m_b+n-1,n);
        
        E0_subspace1_kfe=E0_subspace1_kfe+...
                                sum(dummy_prob(:,index+total_done+1:...
                                           index+total_done+...
                                           number_of_equations),2);
                                           
        total_done=total_done+number_of_equations;
    end
end    

E0_subspace1_pmde=sum(prob(:,total_number_of_equations+1:...
                                    total_number_of_equations+m_a),2);

%error=E0_subspace1_kfe-E0_subspace1_pmde;
%%%%%%%%%%%%%%% Test for E'[N_i^0(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(t,error);
%input('');
%close all

%%%%%%%%%%%%%%% Test for E'[N_i(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1_subspace1_kfe=0;

for l=1:m_a
    index=(l-1)*number_of_equations_per_phase;
    total_done=0;
    
    for n=0:s-1
        number_of_equations=nchoosek(m_b+n-1,n);
        
        E1_subspace1_kfe=E1_subspace1_kfe+n*...
                                sum(dummy_prob(:,index+total_done+1:...
                                           index+total_done+...
                                           number_of_equations),2);
                                           
        total_done=total_done+number_of_equations;
    end
end    

E1_subspace1_pmde=sum(prob(:,total_number_of_equations+m_a+1:...
                       total_number_of_equations+m_a+m_b*m_a),2);

%error=E1_subspace1_kfe-E1_subspace1_pmde;
%%%%%%%%%%%%%%% Test for E'[N_i(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%plot(t,error);
%input('');
%close all


%%%%%%%%%%%%%%% Test for E'[N_i(t)N_j(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E2_subspace1_kfe=0;

for l=1:m_a
    index=(l-1)*number_of_equations_per_phase;
    total_done=0;
    
    for n=0:s-1
        number_of_equations=nchoosek(m_b+n-1,n);
        
        E2_subspace1_kfe=E2_subspace1_kfe+(n^2)*...
                                sum(dummy_prob(:,index+total_done+1:...
                                           index+total_done+...
                                           number_of_equations),2);
                                           
        total_done=total_done+number_of_equations;
    end
end    

E2_subspace1_pmde=sum(prob(:,total_number_of_equations+m_a+m_b*m_a+1:...
                   total_number_of_equations+m_a+m_b*m_a+m_b*m_b*m_a),2);

%error=E2_subspace1_kfe-E2_subspace1_pmde;
%%%%%%%%%%%%%%% Test for E'[N_i(t)N_j(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(t,error);
%input('');
%close all


total_done_subspace1=total_done;
%@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 1 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%



%@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%%%%%%%%%%%%%%%%%%% Test for E'[N^0(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0_subspace2_kfe=0;

for l=1:m_a
    index=(l-1)*number_of_equations_per_phase;
    
    total_done=total_done_subspace1;
    
    for n=s:c
        number_of_equations=nchoosek(m_b+s-1,s);
        
        E0_subspace2_kfe=E0_subspace2_kfe+...
                                sum(dummy_prob(:,index+total_done+1:...
                                           index+total_done+...
                                           number_of_equations),2);
                                           
        total_done=total_done+number_of_equations;
    end
end    

E0_subspace2_pmde=sum(prob(:,total_number_of_equations+m_a+m_b*m_a+...
                               m_b*m_b*m_a+1:total_number_of_equations+...
                               2*m_a+m_b*m_a+m_b*m_b*m_a),2);

%error=E0_subspace2_kfe-E0_subspace2_pmde;
%%%%%%%%%%%%%%%%%%% Test for E'[N^0(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(t,error);
%input('');
%close all



%%%%%%%%%%%%%%%%%%% Test for E'[N(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1_subspace2_kfe=0;

for l=1:m_a
    index=(l-1)*number_of_equations_per_phase;
    
    total_done=total_done_subspace1;
    
    for n=s:c
        number_of_equations=nchoosek(m_b+s-1,s);
        
        E1_subspace2_kfe=E1_subspace2_kfe+n*...
                                sum(dummy_prob(:,index+total_done+1:...
                                           index+total_done+...
                                           number_of_equations),2);
                                           
        total_done=total_done+number_of_equations;
    end
end    

%plot(t,E1_subspace2_kfe);
%input('')
%close all

E1_subspace2_pmde=sum(prob(:,total_number_of_equations+2*m_a+2*m_b*m_a+...
                               m_b*m_b*m_a+1:total_number_of_equations...
                               +3*m_a+2*m_b*m_a+m_b*m_b*m_a),2);

%plot(t,test_variable_pmde);
%input('')
%close all
                           
                           
%error=E1_subspace2_kfe-E1_subspace2_pmde;
%%%%%%%%%%%%%%%%%%% Test for E'[N(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(t,error);
%input('');
%close all

%%%%%%%%%%%%%%%%%%% Test for E'[N^2(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E2_subspace2_kfe=0;

for l=1:m_a
    index=(l-1)*number_of_equations_per_phase;
    
    total_done=total_done_subspace1;
    
    for n=s:c
        number_of_equations=nchoosek(m_b+s-1,s);
        
        E2_subspace2_kfe=E2_subspace2_kfe+(n^2)*...
                                sum(dummy_prob(:,index+total_done+1:...
                                           index+total_done+...
                                           number_of_equations),2);
                                           
        total_done=total_done+number_of_equations;
    end
end    

E2_subspace2_pmde=sum(prob(:,total_number_of_equations+3*m_a+...
                              3*m_b*m_a+2*m_b*m_b*m_a+1:...
                              total_number_of_equations+4*m_a+...
                              3*m_b*m_a+2*m_b*m_b*m_a),2);

error=E2_subspace2_kfe-E2_subspace2_pmde;
%%%%%%%%%%%%%%%%%%% Test for E'[N^2(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot(t,E2_subspace2_pmde);
input('');
close all

E0_kfe=E0_subspace1_kfe+E0_subspace2_kfe;
E1_kfe=E1_subspace1_kfe+E1_subspace2_kfe;
E2_kfe=E2_subspace1_kfe+E2_subspace2_kfe;

E0_pmde=E0_subspace1_pmde+E0_subspace2_pmde;
E1_pmde=E1_subspace1_pmde+E1_subspace2_pmde;
E2_pmde=E2_subspace1_pmde+E2_subspace2_pmde;

E0_error=E0_kfe-E0_pmde;
E1_error=E1_kfe-E1_pmde;
E2_error=E2_kfe-E2_pmde;

%plot(t,E2_error);
%input('')
%close all
%plot(t,E1_pmde);
%input('')
%close all
%plot(t,E2_error);
%input('')
%close all


end