%%%%%%%%%%%%%%%%%%%% PhPhsc Function File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dpdt] = PhPhsc_function(t,prob)
t
global m_a m_b c s total_number_of_equations number_of_equations_per_phase;

[k1,k2,alpha,beta,A1,A2,B1,B2,lambda,mu]=PhPhsc_qparm(t);

%dummy_prob=prob(1:total_number_of_equations,1);

%sum_dummy_prob=sum(dummy_prob);
%dummy_prob=dummy_prob/sum_dummy_prob;

%prob(1:total_number_of_equations,1)=dummy_prob(:,1);

dpdt=zeros(total_number_of_equations+4*m_a+3*m_b*m_a+2*m_b*m_b*m_a,1);

%%%%%%%%%%%%%%%%%%%%%% Test Section (Boundary Moments) %%%%%%%%%%%%%%%%%%%%
%%%%%%%% Subspace 1 terms %%%%%%%%%
E0_bm_n_s_1=zeros(1,m_a);
E1_bm_ni_s_1=zeros(m_b,m_a);
E2_bm_ninj_s_1=zeros(m_b,m_b,m_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dummy_E0_bm_n_s_1=zeros(1,m_a);
%%%%%%%% Subspace 2 terms %%%%%%%%%
E0_bm_n_s=zeros(1,m_a);
E1_bm_ni_s=zeros(m_b,m_a);
E2_bm_ninj_s=zeros(m_b,m_b,m_a);
E3_bm_ninjnk_s=zeros(m_b,m_b,m_b,m_a);
E0_bm_n_c=zeros(1,m_a);
E1_bm_ni_c=zeros(m_b,m_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Those in subspace 1 %%%
E0_n_subspace1=zeros(1,m_a);
E1_ni_subspace1=zeros(m_b,m_a);
E2_ninj_subspace1=zeros(m_b,m_b,m_a);
dummy_E1_n_subspace1=zeros(1,m_a);
dummy_E2_n_subspace1=zeros(1,m_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Those in subspace 2 %%%
E0_n_subspace2=zeros(1,m_a);
E1_ni_subspace2=zeros(m_b,m_a);
E1_n_subspace2=zeros(1,m_a);
E2_ninj_subspace2=zeros(m_b,m_b,m_a);
E2_nin_subspace2=zeros(m_b,m_a);
E2_n_subspace2=zeros(1,m_a);
%dummy_E1_ni_subspace2=zeros(m_b,m_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Kolmogorov Forward Equations %%%%%%%%%%%%%%%%%%%%%%%%%%


for l=1:m_a
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    total_done=0;

    %%%%%%%%%%%%%%%%%%%% P(n_1,...,n_m_b);l,0(t) %%%%%%%%%%%%%%%%%%%%%%%%%%
    container_content=zeros(1,m_b);
    
    for n=0:s-1
          
        container_content(1,1)=n;
        number_of_possibilities=nchoosek(m_b+n-1,n);
        
        for number=1:number_of_possibilities
                        
            current_index=(l-1)*number_of_equations_per_phase+...
                           total_done+number;
            
            %%%%%%%%%
            dpdt(current_index)=-lambda(1,l)*prob(current_index);
            %fine%
            %%%%%%%%%%%
            for j=1:k2
                
                dpdt(current_index)=dpdt(current_index)-...
                                    container_content(1,j)*mu(1,j)*...
                                    (1-B1(j,j))*prob(current_index);               
            end           
            %fine%
            %%%%%%%%%%%%
            for i=1:k1
                influx_index=(i-1)*number_of_equations_per_phase+...
                                                       total_done+number;
                dpdt(current_index)=dpdt(current_index)+...
                                    A1(i,l)*lambda(1,i)*...
                                    prob(influx_index);
            end
            %fine%
            %%%%%%%%%%%%
            if n>0
            for i=1:m_b
                if container_content(1,i)>0
                container_content(1,i)=container_content(1,i)-1;            
                for j=1:k2
                    if j~=i
                            container_content(1,j)=...
                                        container_content(1,j)+1;
                                    
                            influx_index=(l-1)*...
                                    number_of_equations_per_phase+...
                                    total_done;
                            
                            count=0;
                            no_of_containers=m_b;
                            container_no=1;
                            
                            while count~=n
                                  count=count+...
                                        container_content(1,container_no);
                                  balls=n-(count+1);
                                  if balls>=0
                                  influx_index=influx_index+...
                                      nchoosek(no_of_containers+balls-1,...
                                               balls);
                                  else
                                      influx_index=influx_index+1;
                                  end
                                  no_of_containers=no_of_containers-1;
                                  container_no=container_no+1;
                            end
                            
                            dpdt(current_index)=dpdt(current_index)+...
                                B1(j,i)*container_content(1,j)*...
                                mu(1,j)*prob(influx_index);
                            
                            container_content(1,j)=...
                                       container_content(1,j)-1;
                            
                    end     
                end
                container_content(1,i)=container_content(1,i)+1; 
                end 
            end
            end
            %fine%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if n>0
                influx=0;
            for i=1:k1
                
                influx_dummy=0;          
                          
                for j=1:m_b
                    if container_content(1,j)>0
                       container_content(1,j)=container_content(1,j)-1;
                       
                        influx_index=(i-1)*...
                              number_of_equations_per_phase+...
                              total_done-nchoosek(m_b+(n-1)-1,(n-1));                       
                       
                        count=0;
                        no_of_containers=m_b;
                        container_no=1;
                            
                        while count~=(n-1)
                              count=count+...
                                    container_content(1,container_no);
                              balls=(n-1)-(count+1);
                              if balls>=0
                              influx_index=influx_index+...
                              nchoosek(no_of_containers+balls-1,...
                                         balls);
                              else
                                  influx_index=influx_index+1;
                              end
                              no_of_containers=no_of_containers-1;
                              container_no=container_no+1;
                        end
                        
                        if count==0
                            influx_index=influx_index+1;
                        end
                        
                        influx_dummy=influx_dummy+...
                                         beta(1,j)*prob(influx_index);
                        container_content(1,j)=container_content(1,j)+1;
                    end    
                end
                
                influx=influx+A2(i,1)*lambda(1,i)*influx_dummy;
                
            end
            
            dpdt(current_index)=dpdt(current_index)+influx*alpha(1,l);
            
            end
            %fine
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j=1:k2
                container_content(1,j)=container_content(1,j)+1;
                
                influx_index=(l-1)*number_of_equations_per_phase+...
                               total_done+nchoosek(m_b+n-1,n);
                  
                count=0;
                no_of_containers=m_b;
                container_no=1;
                            
                  while count~=n+1
                      count=count+container_content(1,container_no);
                      balls=(n+1)-(count+1);
                      if balls>=0
                          influx_index=influx_index+...
                          nchoosek(no_of_containers+balls-1,balls);
                      else
                          influx_index=influx_index+1;
                      end
                      no_of_containers=no_of_containers-1;
                      container_no=container_no+1;
                  end
                  
                  dpdt(current_index)=dpdt(current_index)+...
                                      B2(j,1)*container_content(1,j)*...
                                      mu(1,j)*prob(influx_index);
                                  
                  container_content(1,j)=container_content(1,j)-1;
            end
            %fine%
            
            %%%%%%%%%% Test Section (To be approximated later) %%%%%%%%%%%%
            dummy_E1_n_subspace1(1,l)=dummy_E1_n_subspace1(1,l)+...
                                      n*prob(current_index);
            
            dummy_E2_n_subspace1(1,l)=dummy_E2_n_subspace1(1,l)+...
                                      (n^2)*prob(current_index);
                                                        
            if n==(s-1)
             
                %dummy_E0_bm_n_s_1(1,l)=dummy_E0_bm_n_s_1(1,l)+...
                %                       prob(current_index);
                
                
                for i=1:m_b
                    E1_bm_ni_s_1(i,l)=E1_bm_ni_s_1(i,l)+...
                                       container_content(1,i)*...
                                       prob(current_index);
                    
                    for j=1:m_b
                        E2_bm_ninj_s_1(i,j,l)=E2_bm_ninj_s_1(i,j,l)+...
                                                container_content(1,i)*...
                                                container_content(1,j)*...
                                                prob(current_index);
                    end
                end
                 
            end
            %%%%%%%%%% Test Section (To be approximated later) %%%%%%%%%%%%
            %%%%%%%%%%%%%%
            if number~=number_of_possibilities
                
            last_non_empty_container=find(container_content,1,'last');
            
            if last_non_empty_container~=m_b
                
            container_content(1,last_non_empty_container)=...    
                container_content(1,last_non_empty_container)-1;
            container_content(1,last_non_empty_container+1)=...
                container_content(1,last_non_empty_container+1)+1;
            
            else
                
            dummy=find(container_content,2,'last');
            second_last_non_empty_container=dummy(1,1);
            balls=container_content(1,last_non_empty_container)+1;
            container_content(1,last_non_empty_container)=0;
            container_content(1,second_last_non_empty_container)=...
                container_content(1,second_last_non_empty_container)-1;
            container_content(1,second_last_non_empty_container+1)=balls;
            
            end
            
            end
            
                
        end
    
        container_content=zeros(1,m_b);
        total_done=total_done+number_of_possibilities;
        
    end
    %fine
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%% P(n_1,...,n_m_b);l,q(t) %%%%%%%%%%%%%%%%%%%%%%%%%
     
    for in_queue=0:c-s
        
        container_content(1,1)=s;
        number_of_possibilities=nchoosek(m_b+s-1,s);
        
        for number=1:number_of_possibilities
        
            current_index=(l-1)*number_of_equations_per_phase+...
                             total_done+number;
                         
            %%%%%%%%%%%%%%
            dpdt(current_index)=-lambda(1,l)*prob(current_index);
            %fine%
            %%%%%%%%%%%%%%
            for j=1:k2
                
                dpdt(current_index)=dpdt(current_index)-...
                                    container_content(1,j)*mu(1,j)*...
                                    (1-B1(j,j))*prob(current_index);
            end
            %fine%
            %%%%%%%%%%%%%%
            for i=1:k1
                
                influx_index=(i-1)*number_of_equations_per_phase+...
                                   total_done+number;
                
                dpdt(current_index)=dpdt(current_index)+...
                                    A1(i,l)*lambda(1,i)*...
                                    prob(influx_index);                 
            end
            %fine%        
            %%%%%%%%%%%%%%
            for i=1:m_b
                if container_content(1,i)>0
                    container_content(1,i)=container_content(1,i)-1;
                for j=1:k2
                    if j~=i
                        if container_content(1,j)<s
                            container_content(1,j)=...
                                             container_content(1,j)+1;
                            
                            influx_index=(l-1)*...
                                number_of_equations_per_phase+total_done;               
                            
                            count=0;
                            no_of_containers=m_b;
                            container_no=1;
                            while count~=s
                                  count=count+...
                                        container_content(1,container_no);
                                  balls=s-(count+1);
                                  if balls>=0
                                  influx_index=influx_index+...
                                      nchoosek(no_of_containers+balls-1,...
                                               balls);
                                  else
                                      influx_index=influx_index+1;
                                  end
                                  no_of_containers=no_of_containers-1;
                                  container_no=container_no+1;
                            end
                            
                            dpdt(current_index)=dpdt(current_index)+...
                                B1(j,i)*container_content(1,j)*...
                                mu(1,j)*prob(influx_index); 
                            
                            container_content(1,j)=...
                                               container_content(1,j)-1;
                        end
                    end
                end
                container_content(1,i)=container_content(1,i)+1;
                end
            end
            %fine%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            influx=0;
            for i=1:k1
                influx_dummy=0;
                
                if in_queue==c-s
                    influx_index=(i-1)*number_of_equations_per_phase+...
                                 total_done+number;
                             
                    influx_dummy=influx_dummy+prob(influx_index);         
                end
                
                if in_queue>0
                    influx_index=(i-1)*number_of_equations_per_phase+...
                                 total_done-nchoosek(m_b+s-1,s)+number;
                    
                    influx_dummy=influx_dummy+prob(influx_index);            
                end
                
                if in_queue==0
                    
                    for j=1:m_b
                        if container_content(1,j)>0
                            container_content(1,j)=...
                                          container_content(1,j)-1;
                                      
                            
                            influx_index=(i-1)*...
                            number_of_equations_per_phase+total_done-...
                            nchoosek(m_b+(s-1)-1,(s-1));
                        
                            count=0;
                            no_of_containers=m_b;
                            container_no=1;
                            
                            while count~=(s-1)
                                count=count+...
                                    container_content(1,container_no);
                                balls=(s-1)-(count+1);
                                if balls>=0
                                 influx_index=influx_index+...
                                 nchoosek(no_of_containers+balls-1,balls);
                                else
                                 influx_index=influx_index+1;
                                end
                                container_no=container_no+1;
                                no_of_containers=no_of_containers-1;
                            end
                            
                            if count==0
                                influx_index=influx_index+1;
                            end
                            
                            influx_dummy=influx_dummy+beta(1,j)*...
                                         prob(influx_index);
                       
                           container_content(1,j)=...
                                          container_content(1,j)+1;
                        end
                    end  
                end
                
                influx=influx+A2(i,1)*lambda(1,i)*influx_dummy;
            end
            
            dpdt(current_index)=dpdt(current_index)+influx*alpha(1,l);
            
            %fine%
           %%%%%%%%%%%%%%%%
           if in_queue<c-s             
               for i=1:k2
                   if container_content(1,i)<s
                       container_content(1,i)=container_content(1,i)+1;
                   for j=1:m_b
                       if j~=i
                           if container_content(1,j)>0
                                
                               container_content(1,j)=...
                                               container_content(1,j)-1;
                               
                               influx_index=(l-1)*...
                               number_of_equations_per_phase+total_done+...
                                             nchoosek(m_b+s-1,s);
                                        
                               count=0;
                               no_of_containers=m_b;
                               container_no=1;
                               
                               while count~=s
                                   count=count+...
                                       container_content(1,container_no);
                                   balls=s-(count+1);
                                   
                                   if balls>=0
                                       influx_index=influx_index+...
                                    nchoosek(no_of_containers+balls-1,...
                                                 balls);
                                   else
                                       influx_index=influx_index+1;
                                   end
                                   no_of_containers=no_of_containers-1;
                                   container_no=container_no+1;  
                               end
                               dpdt(current_index)=dpdt(current_index)+...
                                   B2(i,1)*beta(1,j)*...
                                   container_content(1,i)*mu(1,i)*...
                                   prob(influx_index);
                               
                               container_content(1,j)=...
                                               container_content(1,j)+1;
                                           
                               
                           end
                       end
                   end
                   container_content(1,i)=container_content(1,i)-1;
                   end
               end
                               
               for i=1:k2
                   influx_index=(l-1)*number_of_equations_per_phase+...
                                total_done+nchoosek(m_b+s-1,s)+...
                                number;
                            
                            
                   dpdt(current_index)=dpdt(current_index)+...
                                  B2(i,1)*container_content(1,i)*...
                                  mu(1,i)*beta(1,i)*prob(influx_index);
               end
           end
           %%%%%%%%%%%%%%%%%%%%
           
           %%%%%%%%%%%% Test Section (To be Approximated later) %%%%%%%%%%%           
           if in_queue==0
               
               for i=1:m_b
                   E1_bm_ni_s(i,l)=E1_bm_ni_s(i,l)+...
                                    container_content(1,i)*...
                                    prob(current_index);
                                
                   for j=1:m_b
                       E2_bm_ninj_s(i,j,l)=E2_bm_ninj_s(i,j,l)+...
                                             container_content(1,i)*...
                                             container_content(1,j)*...
                                             prob(current_index);
                                         
                        for k=1:m_b
                            E3_bm_ninjnk_s(j,k,i,l)=...
                                           E3_bm_ninjnk_s(j,k,i,l)+...
                                           container_content(1,i)*...
                                           container_content(1,j)*...
                                           container_content(1,k)*...
                                           prob(current_index);
                        end
                   end
               end
           end
           
           if in_queue==c-s
               
             %  E0_bm_n_c(1,l)=E0_bm_n_c(1,l)+prob(current_index);
               
               for i=1:m_b
                   E1_bm_ni_c(i,l)=E1_bm_ni_c(i,l)+...
                                    container_content(1,i)*...
                                    prob(current_index);
               end
           end
           %%%%%%%%%%%% Test Section (To be Approximated later) %%%%%%%%%%%
                      
           if number~=number_of_possibilities
            last_non_empty_container=find(container_content,1,'last');
            if last_non_empty_container~=m_b
            container_content(1,last_non_empty_container)=...    
                container_content(1,last_non_empty_container)-1;
            container_content(1,last_non_empty_container+1)=...
                container_content(1,last_non_empty_container+1)+1;
            else
            dummy=find(container_content,2,'last');
            second_last_non_empty_container=dummy(1,1);
            balls=container_content(1,last_non_empty_container)+1;
            container_content(1,last_non_empty_container)=0;
            container_content(1,second_last_non_empty_container)=...
                container_content(1,second_last_non_empty_container)-1;
            container_content(1,second_last_non_empty_container+1)=balls;
            end
           end
        end
        
        container_content=zeros(1,m_b);
        total_done=total_done+number_of_possibilities;
    end
    
    
end

%****************** Kolmogorov Forward Equations *************************%

%%%%%%%%%%%%%%%%%%%%%%% Current PMDE values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=1:m_a
    E0_n_subspace1_index=total_number_of_equations + l;
    
    E0_n_subspace2_index=total_number_of_equations + m_a + m_b*m_a +...
                         m_b*m_b*m_a + l;
    
    E1_n_subspace2_index=total_number_of_equations + 2*m_a + 2*m_b*m_a +...
                         m_b*m_b*m_a + l;
                     
    E2_n_subspace2_index=total_number_of_equations+ 3*m_a + 3*m_b*m_a +...
                         2*m_b*m_b*m_a + l;
    
    E0_n_subspace1(1,l)=prob(E0_n_subspace1_index);
    E0_n_subspace2(1,l)=prob(E0_n_subspace2_index);
    E1_n_subspace2(1,l)=prob(E1_n_subspace2_index);
    E2_n_subspace2(1,l)=prob(E2_n_subspace2_index);
    
    for i=1:m_b
        E1_ni_subspace1_index=total_number_of_equations + m_a +...
                              (i-1)*m_a + l;
        
        E1_ni_subspace2_index=total_number_of_equations + 2*m_a + ...
                              m_b*m_a + m_b*m_b*m_a + (i-1)*m_a + l;
                          
        E2_nin_subspace2_index=total_number_of_equations+ 3*m_a +...
                               2*m_b*m_a + 2*m_b*m_b*m_a + (i-1)*m_a + l;
        
        E1_ni_subspace1(i,l)=prob(E1_ni_subspace1_index);
        E1_ni_subspace2(i,l)=prob(E1_ni_subspace2_index);
        E2_nin_subspace2(i,l)=prob(E2_nin_subspace2_index);
        
        for j=1:m_b
            E2_ninj_subspace1_index=total_number_of_equations + m_a +...
                                    m_b*m_a + (l-1)*m_b*m_b + ...
                                    (i-1)*m_b + j;
                                
            E2_ninj_subspace2_index=total_number_of_equations+ 3*m_a +...
                                    2*m_b*m_a + m_b*m_b*m_a + ...
                                    (l-1)*m_b*m_b + (i-1)*m_b + j;
                                
            E2_ninj_subspace1(i,j,l)=prob(E2_ninj_subspace1_index);
            E2_ninj_subspace2(i,j,l)=prob(E2_ninj_subspace2_index);
            
        end
        
    end


end


%%%%%%%%%%%%%%%%%%%%%%% Current PMDE values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%@@@@@@@@@@@@@@@@@@@@@@@@@ Approximations @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%%%%%%%%%%%%%%%%%% Polya Eggenberger Approximations %%%%%%%%%%%%%%%%%%%%%%%

%Subspace1

for l=1:m_a
    
    EN=dummy_E1_n_subspace1(1,l);
    EN2=dummy_E2_n_subspace1(1,l);
    
    EN=EN/E0_n_subspace1(1,l);
    EN2=EN2/E0_n_subspace1(1,l);
    
    if isnan(EN)~=1 && isnan(EN2)~=1 
        
        den=s-1;
           
        Var=EN2-EN*EN;
        eita=1000.;
        x=EN/den;
        epsilon=10^(-4);
          if EN<epsilon
              p=0.;
          else
              p=x;
          end
        theta=p;

        q=1-p;
        pqmin=min(p,q);
        x=-pqmin/(den-1);
        d=EN2-den*EN;
          if(EN<epsilon) 
              gamma=0.;
          elseif(EN>(den-epsilon))
              gamma=0.;
          elseif(Var<epsilon)
              gamma=-1/den+epsilon;
          elseif(abs(d)<epsilon)
              gamma=eita;
          elseif ((EN*(EN+(1.-p)))-EN2)/d<x
              gamma=x+epsilon; 
          else
              gamma=((EN*(EN+(1.-p)))-EN2)/d;
          end

          
          E0_bm_n_s_1(1,l)=PhPhsc_PE(den,1-theta,gamma);
          if E0_bm_n_s_1(1,l)<0
              E0_bm_n_s_1(1,l)=0;
          end
          
          E0_bm_n_s_1(1,l)=E0_bm_n_s_1(1,l)*E0_n_subspace1(1,l);
          
        %if t>1
        %E0_n_subspace1
        %E1_ni_subspace1
        %EN
        %EN2
        %E0_bm_n_s_1(1,l)
        %dummy_E0_bm_n_s_1(1,l)
        %input('');
        %nd
             
    end
end


%Subspace 2
for l=1:m_a
    
    %if c-s==0
     %   E0_bm_n_s(1,:)=E0_n_subspace2(1,:);
      %  E0_bm_n_c(1,:)=E0_n_subspace2(1,:);
    %else
    
    dummy_EN=E1_n_subspace2(1,l)/E0_n_subspace2(1,l);
    dummy_EN2=E2_n_subspace2(1,l)/E0_n_subspace2(1,l);
    
    if isnan(dummy_EN)~=1 && isnan(dummy_EN2)~=1 
        
        EN=dummy_EN-s;
        EN2=dummy_EN2-2*s*EN-s^2;
        
        if EN>=0 && EN2>=0
        
        den=c-s;
           
        Var=EN2-EN*EN;
        eita=1000.;
        x=EN/den;
        epsilon=10^(-4);
          if EN<epsilon
              p=0.;
          else
              p=x;
          end
        theta=p;

        q=1-p;
        pqmin=min(p,q);
        x=-pqmin/(den-1);
        d=EN2-den*EN;
          if(EN<epsilon) 
              gamma=0.;
          elseif(EN>(den-epsilon))
              gamma=0.;
          elseif(Var<epsilon)
              gamma=-1/den+epsilon;
          elseif(abs(d)<epsilon)
              gamma=eita;
          elseif ((EN*(EN+(1.-p)))-EN2)/d<x
              gamma=x+epsilon; 
          else
              gamma=((EN*(EN+(1.-p)))-EN2)/d;
          end

          
          E0_bm_n_s(1,l)=PhPhsc_PE(den,theta,gamma);
          E0_bm_n_c(1,l)=PhPhsc_PE(den,1-theta,gamma);
          
          if E0_bm_n_s(1,l)<0
              E0_bm_n_s(1,l)=0;
          end
             
          if E0_bm_n_c(1,l)<0
              E0_bm_n_c(1,l)=0;
          end
          
          E0_bm_n_s(1,l)=E0_bm_n_s(1,l)*E0_n_subspace2(1,l);
          E0_bm_n_c(1,l)=E0_bm_n_c(1,l)*E0_n_subspace2(1,l);
          
        end
    %end
    
    end
end





%%%%%%%%%%%%%%%%%% Polya Eggenberger Approximations %%%%%%%%%%%%%%%%%%%%%%%



%@@@@@@@@@@@@@@@@@@@@@@@@@ Approximations @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%@@@@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 1 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%%%%%%%%%%%%%%%%%%%%%%%%% E'[N^0(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=1:m_a
        
    dpdt_index=total_number_of_equations+l;
    
    dpdt(dpdt_index)=-lambda(1,l)*E0_n_subspace1(1,l);
    
    for j=1:k1
        dpdt(dpdt_index)=dpdt(dpdt_index)+A1(j,l)*lambda(1,j)*...
                         E0_n_subspace1(1,j);
    end
    
    positive_flux=0;
    negative_flux=0;
    
    for j=1:k1
        positive_flux=positive_flux+A2(j,1)*lambda(1,j)*...
                      E0_n_subspace1(1,j);
                  
        negative_flux=negative_flux+A2(j,1)*lambda(1,j)*...
                      E0_bm_n_s_1(1,j);
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux-...
                     alpha(1,l)*negative_flux;
                 
    for j=1:k2
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(j,1)*mu(1,j)*...
                         E1_bm_ni_s(j,l);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%% E'[N^0(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correct, no problem
%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for l=1:m_a
        dpdt_index=total_number_of_equations + m_a + (i-1)*m_a + l;
        
        dpdt(dpdt_index)=-(lambda(1,l)+mu(1,i))*E1_ni_subspace1(i,l);
        
        for j=1:k1
            dpdt(dpdt_index)=dpdt(dpdt_index)+A1(j,l)*lambda(1,j)*...
                             E1_ni_subspace1(i,j);
        end
        
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B1(j,i)*mu(1,j)*...
                             E1_ni_subspace1(j,l);
        end
        
        positive_flux=0;
        
        for j=1:k1
            positive_flux=positive_flux+A2(j,1)*lambda(1,j)*...
                          E1_ni_subspace1(i,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
        
        positive_flux=0;
        
        for j=1:k1
            positive_flux=positive_flux+A2(j,1)*lambda(1,j)*beta(1,i)*...
                          E0_n_subspace1(1,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
        
        negative_flux=0;
        
        for j=1:k1
            negative_flux=negative_flux+A2(j,1)*lambda(1,j)*...
                          E1_bm_ni_s_1(i,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)-alpha(1,l)*negative_flux;
        
        negative_flux=0;
        
        for j=1:k1
            negative_flux=negative_flux+A2(j,1)*lambda(1,j)*beta(1,i)*...
                          E0_bm_n_s_1(1,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)-alpha(1,l)*negative_flux;
        
        dpdt(dpdt_index)=dpdt(dpdt_index)-B2(i,1)*mu(1,i)*...
                         E1_bm_ni_s(i,l);
                     
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B2(j,1)*mu(1,j)*...
                             E2_bm_ninj_s(i,j,l);
        end
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Not fine, check
%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i^2(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for l=1:m_a
        dpdt_index=total_number_of_equations+m_a+m_b*m_a+(l-1)*m_b*m_b+...
                   (i-1)*m_b+i;
        
        dpdt(dpdt_index)=-lambda(1,l)*E2_ninj_subspace1(i,i,l);
        
        for j=1:k1
            dpdt(dpdt_index)=dpdt(dpdt_index)+A1(j,l)*lambda(1,j)*...
                             E2_ninj_subspace1(i,i,j);
        end
        
        for j=1:k2
            if j~=i
                dpdt(dpdt_index)=dpdt(dpdt_index)+B1(j,i)*mu(1,j)*...
                                 E1_ni_subspace1(j,l)+B1(j,i)*mu(1,j)*...
                                 2*E2_ninj_subspace1(i,j,l);
            end
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+(1-B1(i,i))*mu(1,i)*...
                         (E1_ni_subspace1(i,l)-2*E2_ninj_subspace1(i,i,l));
                     
        positive_flux=0;
        
        for j=1:k1
            positive_flux=positive_flux+A2(j,1)*lambda(1,j)*...
                          E2_ninj_subspace1(i,i,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
        
        positive_flux1=0;
        positive_flux2=0;
        negative_flux1=0;
        negative_flux2=0;
        
        for j=1:k1
            positive_flux1=positive_flux1+A2(j,1)*lambda(1,j)*...
                          E0_n_subspace1(1,j);
            
            positive_flux2=positive_flux2+A2(j,1)*lambda(1,j)*...
                           2*E1_ni_subspace1(i,j);
                      
            negative_flux1=negative_flux1+A2(j,1)*lambda(1,j)*...
                           E0_bm_n_s_1(1,j);
                      
            negative_flux2=negative_flux2+A2(j,1)*lambda(1,j)*...
                           2*E1_bm_ni_s_1(i,j);
                      
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+...
                         alpha(1,l)*beta(1,i)*positive_flux1+...
                         alpha(1,l)*beta(1,i)*positive_flux2-...
                         alpha(1,l)*beta(1,i)*negative_flux1-...
                         alpha(1,l)*beta(1,i)*negative_flux2;
                     
        negative_flux=0;
        
        for j=1:k1
            negative_flux=negative_flux+A2(j,1)*lambda(1,j)*...
                          E2_bm_ninj_s_1(i,i,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)-alpha(1,l)*negative_flux;
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*...
                         E1_bm_ni_s(i,l)-B2(i,1)*mu(1,i)*...
                         2*E2_bm_ninj_s(i,i,l);
                     
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B2(j,1)*mu(1,j)*...
                             E3_bm_ninjnk_s(i,i,j,l);
        end
      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i^2(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Corrected
%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t)N_j(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for j=1:m_b
        if j~=i
            for l=1:m_a
                dpdt_index=total_number_of_equations+m_a+m_b*m_a+...
                           (l-1)*m_b*m_b+(i-1)*m_b+j;
                
                dpdt(dpdt_index)=-(lambda(1,l)+mu(1,i)+mu(1,j))*...
                                 E2_ninj_subspace1(i,j,l);
                             
                for k=1:k1
                    dpdt(dpdt_index)=dpdt(dpdt_index)+A1(k,l)*...
                                     lambda(1,k)*E2_ninj_subspace1(i,j,k);
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)-B1(j,i)*mu(1,j)*...
                                 E1_ni_subspace1(j,l)-B1(i,j)*mu(1,i)*...
                                 E1_ni_subspace1(i,l);
                             
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)+mu(1,k)*...
                                     (B1(k,i)*E2_ninj_subspace1(k,j,l)+...
                                     B1(k,j)*E2_ninj_subspace1(i,k,l));
                end
                
                positive_flux=0;
                
                for k=1:k1
                    positive_flux=positive_flux+A2(k,1)*lambda(1,k)*...
                                  (E2_ninj_subspace1(i,j,k)+...
                                  beta(1,i)*E1_ni_subspace1(j,k)+...
                                  beta(1,j)*E1_ni_subspace1(i,k));
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
                
                negative_flux=0;
                
                for k=1:k1
                    negative_flux=negative_flux+A2(k,1)*lambda(1,k)*...
                                  (E2_bm_ninj_s_1(i,j,k)+...
                                  beta(1,i)*E1_bm_ni_s_1(j,k)+...
                                  beta(1,j)*E1_bm_ni_s_1(i,k));
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)-alpha(1,l)*negative_flux;
                
                dpdt(dpdt_index)=dpdt(dpdt_index)-(B2(i,1)*mu(1,i)+...
                                 B2(j,1)*mu(1,j))*...
                                 E2_bm_ninj_s(i,j,l);
                             
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)+B2(k,1)*mu(1,k)*...
                                     E3_bm_ninjnk_s(i,j,k,l);
                end
                
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t)N_j(t),l,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%@@@@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 1 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%



%@@@@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N^0(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=1:m_a
    dpdt_index=total_number_of_equations+m_a+m_b*m_a+m_b*m_b*m_a+l;
    
    dpdt(dpdt_index)=-lambda(1,l)*E0_n_subspace2(1,l);
    
    for i=1:k1
        dpdt(dpdt_index)=dpdt(dpdt_index)+A1(i,l)*lambda(1,i)*...
                         E0_n_subspace2(1,i);
    end
    
    positive_flux=0;
    
    for i=1:k1
        positive_flux=positive_flux+A2(i,1)*lambda(1,i)*...
                      (E0_bm_n_s_1(1,i)+E0_n_subspace2(1,i));
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
    
    for i=1:k2
        dpdt(dpdt_index)=dpdt(dpdt_index)-B2(i,1)*mu(1,i)*...
                         E1_bm_ni_s(i,l);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N^0(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fine no correction required
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for l=1:m_a
        dpdt_index=total_number_of_equations+2*m_a+m_b*m_a+m_b*m_b*m_a+...
                   (i-1)*m_a+l;
               
        dpdt(dpdt_index)=-(lambda(1,l)+mu(1,i))*E1_ni_subspace2(i,l);
        
        for j=1:k1
            dpdt(dpdt_index)=dpdt(dpdt_index)+A1(j,l)*lambda(1,j)*...
                             E1_ni_subspace2(i,j);
        end
        
        positive_flux1=0;
        positive_flux2=0;
        positive_flux3=0;
        
        for j=1:k1
            positive_flux1=positive_flux1+A2(j,1)*lambda(1,j)*...
                          E1_ni_subspace2(i,j);
            
            positive_flux2=positive_flux2+A2(j,1)*lambda(1,j)*...
                           E1_bm_ni_s_1(i,j);
                       
            positive_flux3=positive_flux3+A2(j,1)*lambda(1,j)*...
                          beta(1,i)*E0_bm_n_s_1(1,j);
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux1+...
                         alpha(1,l)*positive_flux2+...
                         alpha(1,l)*positive_flux3;
        
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B1(j,i)*mu(1,j)*...
                             E1_ni_subspace2(j,l);
        end
        
        positive_flux=0;
        negative_flux=0;
                
        for j=1:k2
            positive_flux=positive_flux+B2(j,1)*mu(1,j)*...
                           beta(1,i)*E1_ni_subspace2(j,l);
                       
            negative_flux=negative_flux+B2(j,1)*mu(1,j)*...
                           (E2_bm_ninj_s(i,j,l)+beta(1,i)*E1_bm_ni_s(j,l));
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+positive_flux-negative_flux;
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*...
                         E1_bm_ni_s(i,l);
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fine, no correction required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=1:m_a
    dpdt_index=total_number_of_equations+2*m_a+2*m_b*m_a+m_b*m_b*m_a+l;
    
    dpdt(dpdt_index)=-lambda(1,l)*E1_n_subspace2(1,l);
    
    for i=1:k1
        dpdt(dpdt_index)=dpdt(dpdt_index)+A1(i,l)*lambda(1,i)*...
                         E1_n_subspace2(1,i);
    end
    
    positive_flux=0;
    
    for i=1:k1
        positive_flux=positive_flux+A2(i,1)*lambda(1,i)*...
                      s*E0_bm_n_s_1(1,i);
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
    
    positive_flux=0;
    
    for i=1:k1
        positive_flux=positive_flux+A2(i,1)*lambda(1,i)*...
                      E1_n_subspace2(1,i);
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
    
    positive_flux=0;
    
    for i=1:k1
        positive_flux=positive_flux+A2(i,1)*lambda(1,i)*...
                      E0_n_subspace2(1,i);
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
    
    negative_flux=0;
    
    for i=1:k1
        negative_flux=negative_flux+A2(i,1)*lambda(1,i)*...
                      E0_bm_n_c(1,i);
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)-alpha(1,l)*negative_flux;
    
    
    for i=1:k2
        dpdt(dpdt_index)=dpdt(dpdt_index)-B2(i,1)*mu(1,i)*...
                         E1_ni_subspace2(i,l);
    end
    
    for i=1:k2
        dpdt(dpdt_index)=dpdt(dpdt_index)-B2(i,1)*mu(1,i)*(s-1)*...
                         E1_bm_ni_s(i,l);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correction required
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i^2(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for l=1:m_a
        dpdt_index=total_number_of_equations+3*m_a+2*m_b*m_a+...
                   m_b*m_b*m_a+(l-1)*m_b*m_b+(i-1)*m_b+i;
               
        dpdt(dpdt_index)=-(lambda(1,l)+2*mu(1,i))*E2_ninj_subspace2(i,i,l);
        
        for j=1:k1
            dpdt(dpdt_index)=dpdt(dpdt_index)+A1(j,l)*lambda(1,j)*...
                             E2_ninj_subspace2(i,i,j);
        end
        
        positive_flux=0;
        
        for j=1:k1
            positive_flux=positive_flux+A2(j,1)*lambda(1,j)*...
                          (E2_ninj_subspace2(i,i,j)+...
                          E2_bm_ninj_s_1(i,i,j)+beta(1,i)*...
                          E0_bm_n_s_1(1,j)+2*beta(1,i)*...
                          E1_bm_ni_s_1(i,j));
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
        
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B1(j,i)*mu(1,j)*...
                             2*E2_ninj_subspace2(i,j,l);
        end
        
        for j=1:k2
            if j~=i
                dpdt(dpdt_index)=dpdt(dpdt_index)+B1(j,i)*mu(1,j)*...
                                 E1_ni_subspace2(j,l);
            end
        end
        
        for k=1:k2
            if k~=i
                dpdt(dpdt_index)=dpdt(dpdt_index)+B1(i,k)*mu(1,i)*...
                                 E1_ni_subspace2(i,l);
            end
        end
        
        positive_flux=0;
        
        for j=1:k2
            if j~=i
                positive_flux=positive_flux+B2(j,1)*mu(1,j)*...
                              (E1_ni_subspace2(j,l)+...
                              2*E2_ninj_subspace2(i,j,l));
            end
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+beta(1,i)*positive_flux;
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*...
                         ((1-beta(1,i))*E1_ni_subspace2(i,l)+...
                         beta(1,i)*2*E2_ninj_subspace2(i,i,l));
                     
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)-B2(j,1)*mu(1,j)*...
                             E3_bm_ninjnk_s(i,i,j,l);
        end
        
        negative_flux=0;
        
        for j=1:k2
            if j~=i
                negative_flux=negative_flux+B2(j,1)*mu(1,j)*...
                              (E1_bm_ni_s(j,l)+2*E2_bm_ninj_s(i,j,l));
            end
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)-beta(1,i)*negative_flux;
        
        dpdt(dpdt_index)=dpdt(dpdt_index)-B2(i,1)*mu(1,i)*(1-beta(1,i))*...
                         (E1_bm_ni_s(i,l)-2*E2_bm_ninj_s(i,i,l));
                
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i^2(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No Correction Required
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t)N_j(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for j=1:m_b
        if j~=i
            for l=1:m_a
                dpdt_index=total_number_of_equations+3*m_a+2*m_b*m_a+...
                           +m_b*m_b*m_a+(l-1)*m_b*m_b+(i-1)*m_b+j;
                       
                dpdt(dpdt_index)=-(lambda(1,l)+mu(1,i)+mu(1,j))*...
                                  E2_ninj_subspace2(i,j,l);
                              
                for k=1:k1
                    dpdt(dpdt_index)=dpdt(dpdt_index)+A1(k,l)*...
                                     lambda(1,k)*E2_ninj_subspace2(i,j,k);
                end
                
                positive_flux1=0;
                positive_flux2=0;
                
                for k=1:k1
                    positive_flux1=positive_flux1+A2(k,1)*lambda(1,k)*...
                                   E2_ninj_subspace2(i,j,k);
                               
                    positive_flux2=positive_flux2+A2(k,1)*lambda(1,k)*...
                                   (E2_bm_ninj_s_1(i,j,k)+beta(1,i)*...
                                   E1_bm_ni_s_1(j,k)+beta(1,j)*...
                                   E1_bm_ni_s_1(i,k));
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)+...
                                 alpha(1,l)*positive_flux1+...
                                 alpha(1,l)*positive_flux2;
                             
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)+B1(k,i)*mu(1,k)*...
                                     E2_ninj_subspace2(k,j,l);
                end
                
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)+B1(k,j)*mu(1,k)*...
                                     E2_ninj_subspace2(i,k,l);
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)-B1(i,j)*mu(1,i)*...
                                 E1_ni_subspace2(i,l)-B1(j,i)*mu(1,j)*...
                                 E1_ni_subspace2(j,l);
                             
                
                positive_flux1=0;
                positive_flux2=0;
                
                for k=1:k2
                    positive_flux1=positive_flux1+B2(k,1)*mu(1,k)*...
                                   E2_ninj_subspace2(k,j,l);
                               
                    positive_flux2=positive_flux2+B2(k,1)*mu(1,k)*...
                                   E2_ninj_subspace2(i,k,l);
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)+beta(1,i)*...
                                 positive_flux1+beta(1,j)*...
                                 positive_flux2;
                             
                dpdt(dpdt_index)=dpdt(dpdt_index)-B2(i,1)*mu(1,i)*...
                                 beta(1,j)*E1_ni_subspace2(i,l)-B2(j,1)*...
                                 mu(1,j)*beta(1,i)*E1_ni_subspace2(j,l);
                             
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)-B2(k,1)*mu(1,k)*...
                                     E3_bm_ninjnk_s(i,j,k,l);
                end
                
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)-B2(k,1)*mu(1,k)*...
                                     beta(1,i)*E2_bm_ninj_s(k,j,l);
                end
                
                for k=1:k2
                    dpdt(dpdt_index)=dpdt(dpdt_index)-B2(k,1)*mu(1,k)*...
                                     beta(1,j)*E2_bm_ninj_s(i,k,l);
                end
                
                dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*...
                                 (E2_bm_ninj_s(i,j,l)+beta(1,j)*...
                                 E1_bm_ni_s(i,l));
                             
                dpdt(dpdt_index)=dpdt(dpdt_index)+B2(j,1)*mu(1,j)*...
                                 (E2_bm_ninj_s(i,j,l)+beta(1,i)*...
                                 E1_bm_ni_s(j,l));
                                 
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t)N_j(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%
%No correction required
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t)N(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m_b
    for l=1:m_a
        dpdt_index=total_number_of_equations+3*m_a+2*m_b*m_a+...
                   2*m_b*m_b*m_a+(i-1)*m_a+l;
               
        dpdt(dpdt_index)=-(lambda(1,l)+mu(1,i))*E2_nin_subspace2(i,l);
        
        for j=1:k1
            dpdt(dpdt_index)=dpdt(dpdt_index)+A1(j,l)*lambda(1,j)*...
                             E2_nin_subspace2(i,j);
        end
        
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B1(j,i)*mu(1,j)*...
                             E2_nin_subspace2(j,l);
        end
        
        positive_flux=0;
        
        for j=1:k1
            positive_flux=positive_flux+A2(j,1)*lambda(1,j)*...
                          (s*E1_bm_ni_s_1(i,j)+...
                          s*beta(1,i)*E0_bm_n_s_1(1,j)+...
                          E2_nin_subspace2(i,j)+E1_ni_subspace2(i,j)-...
                          E1_bm_ni_c(i,j));
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux;
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*...
                         E1_ni_subspace2(i,l);
                     
        for j=1:k2
            dpdt(dpdt_index)=dpdt(dpdt_index)+B2(j,1)*mu(1,j)*...
                             (beta(1,i)*E2_nin_subspace2(j,l)-...
                             E2_ninj_subspace2(i,j,l)-...
                             beta(1,i)*E1_ni_subspace2(j,l)-...
                             beta(1,i)*(s-1)*E1_bm_ni_s(j,l)-...
                             (s-1)*E2_bm_ninj_s(i,j,l));
        end
        
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*(s-1)*...
                         E1_bm_ni_s(i,l);
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N_i(t)N(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correction required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N^2(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=1:m_a
    dpdt_index=total_number_of_equations+3*m_a+3*m_b*m_a+2*m_b*m_b*m_a+l;
    
    dpdt(dpdt_index)=-lambda(1,l)*E2_n_subspace2(1,l);
    
    for i=1:k1
        dpdt(dpdt_index)=dpdt(dpdt_index)+A1(i,l)*lambda(1,i)*...
                         E2_n_subspace2(1,i);
    end
    
    positive_flux1=0;
    positive_flux2=0;
    positive_flux3=0;
    positive_flux4=0;
    negative_flux=0;
    
    for i=1:k1
        positive_flux1=positive_flux1+A2(i,1)*lambda(1,i)*...
                      (s^2)*E0_bm_n_s_1(1,i);
                  
        positive_flux2=positive_flux2+A2(i,1)*lambda(1,i)*...
                       E0_n_subspace2(1,i);
                 
        positive_flux3=positive_flux3+A2(i,1)*lambda(1,i)*...
                       2*E1_n_subspace2(1,i);
                   
        positive_flux4=positive_flux4+A2(i,1)*lambda(1,i)*...
                       E2_n_subspace2(1,i);
                  
        negative_flux=negative_flux+A2(i,1)*lambda(1,i)*(2*c+1)*...
                      E0_bm_n_c(1,i);
    end
    
    dpdt(dpdt_index)=dpdt(dpdt_index)+alpha(1,l)*positive_flux1+...
                     alpha(1,l)*positive_flux2+...
                     alpha(1,l)*positive_flux3+...
                     alpha(1,l)*positive_flux4-...
                     alpha(1,l)*negative_flux;
                 
    for i=1:k2
        dpdt(dpdt_index)=dpdt(dpdt_index)+B2(i,1)*mu(1,i)*...
                         (E1_ni_subspace2(i,l)-2*E2_nin_subspace2(i,l)-...
                         ((s-1)^2)*E1_bm_ni_s(i,l));
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E'[N^2(t),l,2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%@@@@@@@@@@@@@@@@@@@@@@@@@@@ Subspace 2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%

%*************************************************************************%
%@@@@@@@@@@@@@@@@@@@@@@@@ Moment Equations @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%
end