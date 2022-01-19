%%% Filename: ED_SpinHalf_EScolumns_HeisenbergChains.m | Description: ... | Author: Seyed Nariman Saadatmand | Contact: n.saadatmand@griffith.edu.au | Created in: 29/May/2018 | VER="1.01.06"

%%% Notes: 
% - We use a bit-representation of the Sz-computational basis as described in [Sandvik, AIP Conference 
%   Proceedings, 1297(1):135–338, 2010].
% - 1D particle number, i.e. physical locations, is running from 0 to N-1 (the code is restricted to even N).
% - Correspondingly, states numbering, i.e. the "bits", are running from 0 to 2^N-1 (however, note that Matlab 
%   indexing for arrays always starts from 1).
% - In this code, spin rotational SU(2) symmetry (implemented by fixing the total magnetizations in Z-direction)
%   can be preserved as a good quantum number for the Heisenberg model and one can be read as (symmetry charge)  
%   q-values in parallel with the eigenvectors.
% - The lattice has a 1D structure and can be considered with OBC or PBC (in the latter case the translational 
%   invariance can be also used to find states MOMENTA).
% - The code can now perform Lanczos approach (in addition to the full diagonalization scheme) to find a subset of 
%   eigenvalues/vectors for the Hamiltonians having a sparse matrix form.
% - The code will automatically produce an output file containing a list of useful observables (see below).
% - FOR LATER: for the 'classical' Hamiltonians, the code effectively performs NO diagonalization, but just checks 
%   and report the lowest energy states in the Sz-computational basis.   
% - FOR LATER: the default output name is set to LotusCCT-ED_SpinHalf_chain-[time][day][month][year].out. 


function [STATES, U_H, H_diag] = ED_Spin1_InfRange_HeisenbergChains(N, TolE, q, lanczos, n_eig, MomentumStates)


 %%% printing the welcome message to STDOUT:
 fprintf('\nA Matlab function | ''Spin-1 infinite-range XX Heisenberg chain'' project | VER=''1.01.16''\n');
 fprintf('Copyright (C) Seyed Nariman Saadatmand 2018\n');
 fprintf('Contact: n.saadatmand@griffith.edu.au\n');
 fprintf('NOTE1: This program comes with NO SUPPORT; this is a free software, and you are welcome to\n'); 
 fprintf('redistribute it under certain conditions. Please contact me for further details.\n');
 fprintf('NOTE2: Authors of research publications making use of this software should contact me to discuss citation/acknowledgement conditions\n');
 fprintf('prior to the publication of their works.\n');
 fprintf('USAGE: Lotus_CompCollection_ED_SpinHalf_chain(N, TolE, q, BuiltinSym, n_EigReport, BC, MomentumStates, J_Sz, J_Szz, J_SS, J_3N, J_4N, J_LR, alpha)\n');
 fprintf('OPTIONS DESCRIPTION: ...\n');
 fprintf('EXAMPLE RUN: [STATES, U_H, H_diag] = LotusED_SpinHalf_LongRangeHeisenbergChain(6, 1e-12, 0, ''SU(2)'', ''YES'', 20, ''PBC'', ''YES'', 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0)\n');
 fprintf('PREREQUISITES: this program requires the presence of the functions ''tritget'',''QutritCycle'', ''FlipUpDown'', ''FlipDownUp'', and ''FlipS'' in the CWD.\n\n');

 %%% initial values/settings for global usage:
 %DIR='.';
 DIR=strcat(getenv('HOME'),'/Dropbox/AcademiaJobs-eDesktop/MyPapers/InteriorPoint--InfiniteRangeHams/MatlabCollection');
 setenv('EDITOR','vim');
 
 %%% setting up the output files:
 MainOutName = strcat('EnergySpectrum_q',num2str(q),'-EDresults-L',num2str(N),'-InfRangeHeis.out');
 if exist(fullfile(DIR,MainOutName),'file')
   fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,MainOutName));
   FileID_main = fopen( fullfile(DIR,MainOutName) , 'at');
   if FileID_main==-1
     error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,MainOutName));
   end
 else
   %error('ERROR: this file does not exist: %s', fullfile(DIR,filename4));   
   edit(fullfile(DIR,MainOutName));
   FileID_main = fopen( fullfile(DIR,MainOutName) , 'at');
   if FileID_main==-1
     error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,MainOutName));
   end
   fprintf(FileID_main,'#EnergyIndex\t#energy_per_site\t#momentum_per_PI\t#Mz1_per_site\t#Mz2_per_site\n');   % printing the file header  
 end

 
 %STATES = zeros(3^N,1);
 STATES = int64.empty(1,0);   % the cost of constructing such a vector should scale as O(2^N/sqrt(N)) ...
 qq=0;
 for pp = 1 : 3^N
    
   Sz_total=0;  
   for kk = 1 : N        
     Sz_total = Sz_total + ( double(tritget(pp-1,kk)) - 1 );
   end
     
    %if strcmp(BuiltinSym,'NONE') || ( strcmp(BuiltinSym,'U(1)') && (Sz_total==q) )
    if Sz_total==q
      qq = qq+1;  
      STATES(qq) = pp-1;  
    end    
    
 end 
 STATES = STATES';
  
            
 %%% Constructing the Hamiltonian:
 fprintf('NOTE: constructing the Hamiltonian according to the inserted parameters ...\n');
 %fprintf('J_Sz=%.2f as the coefficient of Sz-field term\n', J_Sz);
 %fprintf('J_Szz=%.2f as the coefficient of Sz*Sz term\n', J_Szz);
 %fprintf('J_SS=%.2f as the coefficient of the NN Heisenberg term\n', J_SS);
 %fprintf('J_3N=%.2f as the coefficient of the NNN Heisenberg term\n', J_3N);
 %fprintf('J_4N=%.2f as the coefficient of the NNNN Heisenberg term\n', J_4N);
 %fprintf('J_LR=%.2f as the coefficient of the LR Heisenberg term (supposing PBC)\n', J_LR);
 %fprintf('alpha=%.2f as the power of the LR Heisenberg term (supposing PBC)\n', alpha);

 H = zeros( size(STATES,1) );

 for ss = 1:size(STATES,1)     
     
  for ii = 0 : N-2

   %fprintf('DEBUG: before FIELD and for (ss=%i,ii=%i), we have H(ss,ss)=%1.4f \n', ss, ii, H(ss,ss)); %___DEBUG___%
   %H(ss,ss) = H(ss,ss) + ( double(tritget(STATES(ss),ii+1)) - 1/2 )*J_Sz;
   %fprintf('DEBUG: after FIELD and for (ss=%i,ii=%i), we have H(ss,ss)=%1.4f \n', ss, ii, H(ss,ss)); %___DEBUG___%

   for jj = mod(ii + 1,N) : N-1 

    %if ( ii < N-1 ) || ( (ii==N-1) && strcmp(BC,'PBC') )     

      if ( tritget(STATES(ss),ii+1) == tritget(STATES(ss),jj+1) )   

        %if ( mod(ii + 1,N) == jj ) 
        %  H(ss,ss) = H(ss,ss) + (1/4)*(J_Szz+J_SS);
        %elseif ( mod(ii + 2,N) == jj )
        %  H(ss,ss) = H(ss,ss) + (1/4)*J_3N;
        %elseif ( mod(ii + 3,N) == jj )
        %  H(ss,ss) = H(ss,ss) + (1/4)*J_4N;
        %end
        
        if ( tritget(STATES(ss),ii+1) == 1 ) && ( tritget(STATES(ss),jj+1) == 1 )
          next_index_UpDown = find( STATES == FlipUpDown( N, STATES(ss), ii, jj ) );  
          H(ss,next_index_UpDown) = H(ss,next_index_UpDown) + 2; 
          next_index_DownUp = find( STATES == FlipDownUp( N, STATES(ss), ii, jj ) ); 
          H(ss,next_index_DownUp) = H(ss,next_index_DownUp) + 2; 
        end
        
      else 
        
        if ( tritget(STATES(ss),ii+1) == 1 ) || ( tritget(STATES(ss),jj+1) == 1 ) 
          next_index = find( STATES == FlipS( N, STATES(ss), ii, jj ) );      
          H(ss, next_index) = H(ss, next_index) + 2;
        else
          next_index = find( STATES == Zeroing( N, STATES(ss), ii, jj ) );      
          H(ss, next_index) = H(ss, next_index) + 2;  
        end  
  
      end      
   
    %end

   end    

  end

 end  

 
%%%%%%%%%%%%%%%%%%%% sparse/full diagonalization (if necessary) %%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if strcmp(lanczos,'YES')
  fprintf('NOTE: starting the Lanczos digonalization of the Hamiltonain using ''eigs'' function; this may take a long time ...\n');
  opts.tol = TolE;
  opts.maxit = 500;
  %opts.p = max(20,2*n_eig);
  [U_H,H_diag] = eigs(H,n_eig,'sa',opts);
  %H_diag = rot90(fliplr(H_diag),-1);
  %U_H = fliplr(U_H);
 else
  fprintf('NOTE: starting the full digonalization of the Hamiltonain using ''eig'' function; this may take a long time ...\n');
  [U_H,H_diag] = eig(H);
 end
 H_array = diag(H_diag)';
 
 
 disp(U_H); %DEBUG
 %disp(H); %DEBUG
 
 
 %%% constructing DegE_array:
 i_energy = -1;   % this would be the index for the unique energies in H_array
 E_old = -9999.0;
 DegE_array = zeros(1,size(H_array,2));

 for ii = 1 : size(H_array,2) 
 
    E = H_array(ii);
    
    if ( abs(E-E_old) > TolE )
      i_energy = i_energy+1;  
      E_old = E;
    end

    DegE_array(ii) = i_energy;

 end  

    
%%%%%%%%%% If we have PBC and 'MomentumStates==YES', calculate the lattice momenta for the states %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if strcmp(MomentumStates,'YES')

   fprintf('NOTE: calculating the lattice momenta of the states ...\n');
   
   %%% Finding the 'big' T-operator:

   T = zeros( size(STATES,1) );
 
   for ss = 1 : size(STATES,1)    
      %StateT_index = find( STATES == CycleBits( STATES(ss), N ) );
      StateT_index = ( STATES == CycleBits( STATES(ss), N ) );
      T(ss,StateT_index) = 1;   
   end 
 
   T_expectation = U_H' * T * U_H;

   %%% Finding the T-operator for a block of all-equal eigenvalues in H:

   i_energy = 0;
   TrialEnergy_array = -9999.0*ones(1,size(H_array,2));
   TT = double.empty(1,0);

   for ii = 1 : size(H_array,2)
   
     E = H_array(ii);
  
     if all( abs( E - TrialEnergy_array ) > TolE )
    
       i_energy = i_energy+1;  
       TrialEnergy_array(i_energy) = E;
       jj = 0;
       T_IndexInfo = double.empty(1,0);

       for kk = ii : size(H_array,2)
        E_test = H_array(kk);   
        if ( abs(E_test - E) < TolE )   
          jj = jj+1;
          T_IndexInfo(jj) = kk;
        end  
       end
  
       TT{ii}=T_IndexInfo;

     end
   
   end

   momentum = zeros( size(H_array,2) );
   [~, size_TT] = size(TT);

   %%% diagonalizing the T-operator for each unique block of all-equal energies:

   for pp = 1 : size_TT 
   
      if ~isempty(TT{1,pp})
  
        T_block = double.empty;    
 
        %%% forming the small-block T-matrices:
        for q1 = 1 : size( TT{1,pp}, 2)
         for q2 = 1 : size( TT{1,pp}, 2)   
          T_block(q1,q2) = T_expectation( TT{1,pp}(q1), TT{1,pp}(q2) );
         end 
        end

        [~,T_diag] = eig(T_block);

        %%% now, filling the momentum states:
        for qq = 1 : size(T_diag,1)
           momentum( TT{1,pp}(qq) , TT{1,pp}(qq) ) = T_diag(qq,qq); 
        end 

      end

   end

   p_angle = angle(momentum)/pi;

 else

  fprintf('NOTE: calculation of lattice momenta is NOT possible/requested ...\n');
  p_angle = -9999.0*eye( size(STATES,1) );

 end   
 

 %%% Calculating the magnetization order parameters for the lowest-energy states:
 fprintf('NOTE: Calculating the magnetization order parameters for the lowest-energy states ...\n\n');
 
 Mz1 = zeros( size(STATES,1) );
 Mz2 = zeros( size(STATES,1) );
 %Mz3 = zeros( size(STATES,1) );

 for ss = 1 : size(STATES,1)

   Mz1(ss,ss) = 0.0;
   Mz2(ss,ss) = 0.0;
   %Mz3(ss,ss) = 0.0;  

   for kk = 1 : N        
     Mz1(ss,ss) = Mz1(ss,ss) + ( double(tritget(STATES(ss),kk)) - 1 );
     Mz2(ss,ss) = Mz2(ss,ss) + (-1)^kk*( double(tritget(STATES(ss),kk)) - 1 );   
     %Mz3(ss,ss) = Mz3(ss,ss) + (2*tritget(STATES(ss),kk)-1)/4;
   end

 end 

 Mz1_expectation = (1/N) * U_H' * Mz1 * U_H;
 Mz2_expectation = (1/N) * U_H' * Mz2 * U_H;
 %Mz3_expectation = (1/N) * U_H' * Mz3 * U_H;
 %Mz3_module = mod(Mz3_expectation);
  

 %%% Printing onto the output file ...
 fprintf('NOTE: Printing observables onto the main output file up to ''n_eig'' ...\n');

 for ii = 1 : min(n_eig,size(STATES,1))
    fprintf('%i\t%1.16f\t%1.2f\t%1.16f\t%1.16f\n', DegE_array(ii), H_array(ii)/N, p_angle(ii,ii), Mz1_expectation(ii,ii),  Mz2_expectation(ii,ii)); 
    fprintf(FileID_main,'%i\t%1.16f\t%1.2f\t%1.16f\t%1.16f\n', DegE_array(ii), H_array(ii)/N, p_angle(ii,ii), Mz1_expectation(ii,ii),  Mz2_expectation(ii,ii)); 
 end

 fclose(FileID_main);

end