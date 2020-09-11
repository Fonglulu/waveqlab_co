#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 11:48:03 2019

@author: shilu fred 
"""

import numpy as np





def backward_difference(u_old, num_x, h, order):
    
    u_new = np.zeros(np.shape(u_old))
    
    m = num_x-1
    
    
    if order == 32:
        
        
        u_new[0] = -u_old[0]+ u_old[1]
        u_new[1] = -0.6923076923076923*u_old[0] +0.3846153846153846*u_old[1]+\
                     0.3076923076923077*u_old[2]
        u_new[m-1] = 0.1538461538461539*u_old[m-3]- 0.9230769230769231*u_old[m-2]+\
                      0.3846153846153846*u_old[m-1] + 0.3846153846153846*u_old[m]
        u_new[m] =  0.4*u_old[m-2]-1.8*u_old[m-1]+1.4*u_old[m]

        for step_i in range(2,m-1):

            u_new[step_i] = 1/6 * u_old[step_i-2] - u_old[step_i-1] + 0.5* u_old[step_i] +1/3 * u_old[step_i+1]
            
            
            
            
    if order == 54:
        
        u_new[0] =  -1.410358565737052*u_old[0] + 1.731075697211155*u_old[1] \
                   -0.231075697211155*u_old[2]-0.089641434262948*u_old[3]
                   
        u_new[1] = -0.524526198439242 * u_old[0] + 0.073578595317726  * u_old[1]\
                   + 0.426421404682274* u_old[2] + 0.024526198439242 * u_old[3]
        
        u_new[2] = 0.148499210110585 * u_old[0]  - 0.888625592417062 *u_old[1]\
                   + 0.274881516587678 * u_old[2] + 0.522116903633491*u_old[3]\
                   - 0.056872037914692 * u_old[4]
                   
        u_new[3] = 0.014208389715832* u_old[0] + 0.116373477672530 * u_old[1]\
                   - 0.885656292286874 * u_old[2] + 0.316644113667118 * u_old[3]\
                   + 0.487144790257104 * u_old[4] - 0.048714479025710  *u_old[5]
                   
                   
        u_new[m-3] = -0.032476319350474  *u_old[m-6] + 0.243572395128552  * u_old[m-5]\
                     -0.974289580514208 *u_old[m-4] + 0.316644113667118 * u_old[m-3]\
                     + 0.447225981055480 * u_old[m-2]  + 0.029769959404601 *u_old[m-1]\
                     - 0.030446549391069 * u_old[m]
                     
        u_new[m-2] = -0.037914691943128  * u_old[m-5] + 0.284360189573460 * u_old[m-4]\
                     -1.033965244865719 * u_old[m-3] + 0.274881516587678  * u_old[m-2]\
                     + 0.604265402843602 * u_old[m-1] - 0.091627172195893 * u_old[m]
                     
        u_new[m-1] = -0.026755852842809  * u_old[m-4] + 0.095875139353400* u_old[m-3]\
                     -0.627090301003345  * u_old[m-2] + 0.073578595317726 * u_old[m-1]\
                     + 0.484392419175028 * u_old[m]
        
        u_new[m] = 0.041832669322709 * u_old[m-3] +0.374501992031873 *u_old[m-2]\
                   -1.874501992031873 * u_old[m-1] + 1.458167330677291 * u_old[m]
                   
                   
        for i in range(4, m-3):
        
            u_new[i] = -1/float(30) * u_old[i -3] + 1/float(4)* u_old[i-2] - u_old[i -1]\
                       +1/float(3) * u_old[i] + 1/float(2) *u_old[i+1] - 1/float(20) * u_old[i+2]
                       
                       
                       
    if order  == 6:
        
        u_new[0] = - 1.580183696207410*u_old[0] + 2.039450367044345*u_old[1] \
                   - 0.189297839436616*u_old[2] - 0.366971721882124*u_old[3] \
                   + 0.044953974933765*u_old[4] + 0.052048915548040*u_old[5]
                   
                   
        u_new[1] = - 0.473094335557797*u_old[0] + 0.029551468764759*u_old[1]\
                   + 0.279404147185598*u_old[2] + 0.215422101432619*u_old[3]\
                   - 0.021790841692085*u_old[4] - 0.029492540133095*u_old[5]
                   
                   
        u_new[2] =  0.140255933237354*u_old[0] - 0.885608401570530*u_old[1] \
                   + 0.306540940575249*u_old[2] + 0.491468255323896*u_old[3] \
                   - 0.061405392278187*u_old[4] + 0.008748664712219*u_old[5]
                   
                   
        u_new[3] = 0.076418487423426*u_old[0] - 0.118346128899150*u_old[1] \
                  -0.624133691970994*u_old[2] + 0.318292975073623*u_old[3] \
                  +0.327107204245208*u_old[4] + 0.020661154127888*u_old[5]
                  
        
        u_new[4] = - 0.008331008040049*u_old[0] - 0.047235828122731*u_old[1]\
                   + 0.402190050630647*u_old[2] - 1.303115093827634*u_old[3]\
                   + 0.593528379527062*u_old[4] + 0.399693508760138*u_old[5]\
                   - 0.036730008927433*u_old[6]
                   
                 
        u_new[5] = - 0.016775726857039*u_old[0] + 0.054367610816546*u_old[1]\
                   - 0.120254802392339*u_old[2] + 0.406675342006739*u_old[3]\
                   - 1.257688453771463*u_old[4] + 0.572337434946968*u_old[5]\
                   + 0.394187558455186*u_old[6] - 0.032848963204599*u_old[7]
                   
                   
        u_new[m-5] = 0.016424481602299*u_old[m-9] - 0.131395852818395*u_old[m-8]\
                   + 0.492734448068983*u_old[m-7] - 1.313958528183954*u_old[m-6]\
                   + 0.572337434946968*u_old[m-5] + 0.357460227911152*u_old[m-4]\
                   + 0.025427590054054*u_old[m-3] + 0.005338567444210*u_old[m-2]\
                   - 0.040531470682718*u_old[m-1] + 0.016163101657400*u_old[m]
                   
        u_new[m-4] = 0.018365004463716*u_old[m-8] - 0.146920035709731*u_old[m-7]\
                   + 0.550950133911491*u_old[m-6] - 1.406282074938901*u_old[m-5]\
                   + 0.593528379527062*u_old[m-4] + 0.450132219590077*u_old[m-3]\
                   - 0.041897570295664*u_old[m-2] - 0.033485250803205*u_old[m-1]\
                   + 0.015609194255155*u_old[m]
                   
                  
        u_new[m-3] = 0.013345690454124*u_old[m-7] - 0.106765523632994*u_old[m-6]\
                   + 0.330443502642251*u_old[m-5] - 0.946962507015985*u_old[m-4]\
                   + 0.318292975073623*u_old[m-3] + 0.243684445410829*u_old[m-2]\
                   + 0.240558007429252*u_old[m-1] - 0.092596590361100*u_old[m]
                   
                   
        u_new[m-2] = 0.026915887850467*u_old[m-6] - 0.197069524204597*u_old[m-5]\
                   + 0.589452745232692*u_old[m-4] - 1.258766829227482*u_old[m-3]\
                   + 0.306540940575249*u_old[m-2] + 0.629259919521358*u_old[m-1]\
                   - 0.096333139747688*u_old[m]
                   
                   
        u_new[m-1] = 0.039560344515975*u_old[m-5] - 0.030739159126096*u_old[m-4]\
                   - 0.105980141988699*u_old[m-3] - 0.393228064437077*u_old[m-2]\
                   + 0.029551468764759*u_old[m-1] + 0.460835552271137*u_old[m]
                   
                   
        u_new[m] = - 0.054021710000147*u_old[m-5] - 0.023993033880123*u_old[m-4]\
                   + 0.302855902188632*u_old[m-3] + 0.275607596716317*u_old[m-2]\
                   - 2.093702214477299*u_old[m-1] + 1.593253459452621*u_old[m]
                   
                   
        for i in range(6,m-5):
            
            u_new[i] = 0.016666666666667*u_old[i-4] - 0.133333333333333*u_old[i-3] + 0.5*u_old[i-2]\
                      - 1.333333333333333*u_old[i-1] + 0.583333333333333*u_old[i] + 0.4*u_old[i+1]\
                      - 0.033333333333333*u_old[i+2]
                   
                    
        return u_new/h
    
    
    
                   
                   
                  
                       
                   
                   
                       
    if order == 7:
            
            u_new[0] = -1.581527971597573 * u_old[0] + 2.045543388652316 * u_old[1]  -0.200227171966867 * u_old[2]\
                       - 0.357299100037565 * u_old[3] + 0.040746019354332 * u_old[4] + 0.052764835595357 * u_old[5]
                       
            u_new[1] = -0.468965247789993 * u_old[0] + 0.012686575154571 * u_old[1] + 0.305572843948314* u_old[2]\
                       + 0.196814495127562 * u_old[3] - 0.016267583768385 * u_old[4] - 0.029841082672069 * u_old[5]
                       
                       
            u_new[2] = 0.119950629785382 * u_old[0] - 0.790802931381985 * u_old[1] + 0.130372094340790 * u_old[2]\
                       + 0.654195007415723 * u_old[3] - 0.136047721252784 * u_old[4] + 0.022332921092875 * u_old[5]
                       
            u_new[3] = 0.083677769468176 * u_old[0] - 0.167666797363523 *  u_old[1] - 0.491788539658796 * u_old[2]\
                       + 0.136933409573552 * u_old[3] + 0.461432440844702 * u_old[4] - 0.030243581766320 * u_old[5]\
                     + 0.007655298902209 * u_old[6]
                       
            u_new[4] = - 0.011021324291367 *u_old[0] - 0.009614356523787 * u_old[1] + 0.247302955053649 * u_old[2]\
                       - 0.995508274694755 * u_old[3] + 0.253703695582764 * u_old[4]  + 0.614572168859249 * u_old[5]\
                       - 0.109901691773728 * u_old[6] + 0.010466827787974 * u_old[7]
                       
            u_new[5] = -0.016680563512292 * u_old[0] + 0.046892232935029 * u_old[1] - 0.063012072625329 * u_old[2]\
                       + 0.222370039269531 * u_old[3] - 0.937327875498694 * u_old[4] + 0.245402632435032 * u_old[5]\
                       + 0.591558939080254 * u_old[6]  - 0.098593156513376 * u_old[7] + 0.009389824429845 * u_old[8]
                       
            u_new[m-5] = 0.007042368322384 * u_old[m-9] - 0.065728771008917 * u_old[m-8] + 0.295779469540127 * u_old[m-7]\
                        -  0.985931565133756 * u_old[m-6] + 0.245402632435032 * u_old[m-5] + 0.551334643308827 * u_old[m-4]\
                        - 0.037096124729170 * u_old[m-3]  + 0.013739135291768 * u_old[m-2] - 0.040959674614960 * u_old[m-1]\
                        + 0.016417886588667 *u_old[m]
                        
            u_new[m-4] = 0.007850120840981 * u_old[m-8] - 0.073267794515818 * u_old[m-7] + 0.329705075321183 * u_old[m-6]\
                        - 1.044838434095626 * u_old[m-5] + 0.253703695582764 * u_old[m-4] + 0.630900759826953 * u_old[m-3]\
                        - 0.093295940554189 * u_old[m-2] - 0.024889867269617 * u_old[m-1] + 0.014132384863370 * u_old[m]
                        
            u_new[m-3] = 0.005741474176657 * u_old[m-7] - 0.053587092315462 * u_old[m-6] + 0.181292965616422 * u_old[m-5]\
                         - 0.728101537236212 * u_old[m-4] + 0.136933409573552 * u_old[m-3] + 0.328114765554538 * u_old[m-2]\
                         + 0.220243795689184 * u_old[m-1] - 0.090637781058677 * u_old[m]
                         
            u_new[m-2] =  0.011447347501192 * u_old[m-6] - 0.102425925355228 * u_old[m-5] + 0.360626660648614 * u_old[m-4]\
                          - 0.980527672399361 * u_old[m-3] + 0.130372094340790 * u_old[m-2] + 0.681777727691031 * u_old[m-1]\
                         - 0.101270232427039 *u_old[m]
                         
            u_new[m-1] = 0.034163235251410 * u_old[m-5] - 0.006283775981431 * u_old[m-4] - 0.149830581921711 * u_old[m-3]\
                        - 0.354437950860383 * u_old[m-2] + 0.012686575154571 * u_old[m-1] + 0.463702498357543 * u_old[m]
                        
                        
            u_new[m] = - 0.053609043198748 * u_old[m-5] - 0.031776313568303 * u_old[m-4] + 0.329862352927359 * u_old[m-3]\
                        + 0.237161254615221 * u_old[m-2] - 2.068759097745568 * u_old[m-1] + 1.587120846970038 * u_old[m]
                        
                        
            for i in range(6, m-5):
                
                u_new[i] =  0.007142857142857 * u_old[i-4] - 0.066666666666667 * u_old[i-3] + 0.3 * u_old[i-2]\
                            - u_old[i-1] + 0.25 *u_old[i] + 0.6 * u_old[i+1] - 0.1 * u_old[i+2] + 0.009523809523810 * u_old[i+3]
                
                
                      
                    
                  
    return u_new/h
                  
        



def forward_difference(u_old, num_x, h, order):
    
    u_new = np.zeros(np.shape(u_old))
    
    m = num_x-1
    
    if order == 32:
        
  
        
        u_new[0] = -1.4*u_old[0]+1.8*u_old[1]-0.4*u_old[2]
        u_new[1] = -0.3846153846153846*u_old[0]- 0.3846153846153846*u_old[1] +\
                    0.9230769230769231*u_old[2]- 0.1538461538461539*u_old[3]
        
        u_new[m-1] = -0.3076923076923077*u_old[m-2] - \
                     0.3846153846153846*u_old[m-1] + 0.6923076923076923*u_old[m]
        u_new[m] = -u_old[m-1]+u_old[m]

        for step_i in range(2,m-1):

            u_new[step_i] = -1/3 * u_old[step_i-1] - 0.5*u_old[step_i] + u_old[step_i+1] -1/6 *u_old[step_i+2]
            
            
            
    
    if order == 54:
        
        u_new[0] =  -1.458167330677291 * u_old[0] + 1.874501992031873 * u_old[1] - 0.374501992031873 *u_old[2]\
                    - 0.041832669322709 * u_old[3]
                    
        u_new[1] =  - 0.484392419175028 * u_old[0] - 0.073578595317726  *u_old[1] + 0.627090301003345 *u_old[2]\
                    - 0.095875139353400 * u_old[3] + 0.026755852842809 * u_old[4]
                    
        u_new[2] = 0.091627172195893  * u_old[0] - 0.604265402843602  * u_old[1] - 0.274881516587678  * u_old[2]\
                   + 1.033965244865719 * u_old[3] - 0.284360189573460 * u_old[4] + 0.037914691943128  * u_old[5]\
            
            
        u_new[3] = 0.030446549391069 * u_old[0] - 0.029769959404601 *u_old[1] - 0.447225981055480 *u_old[2]\
                  - 0.316644113667118 *u_old[3] +  0.974289580514208  * u_old[4] - 0.243572395128552 * u_old[5]\
                  + 0.032476319350474 * u_old[6]
                  
                  
                  
                  
        u_new[m-3] = 0.048714479025710 * u_old[m-5]  - 0.487144790257104  *u_old[m-4] - 0.316644113667118 *u_old[m-3]\
                      + 0.885656292286874 * u_old[m-2] - 0.1163734776725301 * u_old[m-1]  - 0.014208389715832 *u_old[m]
                      
        
        
        u_new[m-2] = 0.056872037914692 * u_old[m-4] - 0.522116903633491 * u_old[m-3] - 0.274881516587678 * u_old[m-2]\
                     + 0.888625592417062* u_old[m-1] - 0.148499210110585 * u_old[m]
                     
                     
        u_new[m-1] = - 0.024526198439242 * u_old[m-3] - 0.426421404682274 * u_old[m-2] - 0.073578595317726 *u_old[m-1]\
                     + 0.524526198439242 * u_old[m]
                     
                     
        u_new[m] = 0.089641434262948 * u_old[m-3] + 0.231075697211155* u_old[m-2] - 1.731075697211155 *u_old[m-1]\
                   +1.410358565737052 * u_old[m]


        for i in range(4, m-3):
        
            u_new[i] = 1/float(20) * u_old[i -2] - 1/float(2) * u_old[i-1] - 1/float(3) * u_old[i]\
                           + u_old[i+1 ] - 1/float(4) *u_old[i+2] + 1/float(30) * u_old[i+3]
                        
                           
                        
    if order == 6 :
        
        u_new[0] = - 1.593253459452621*u_old[0] + 2.093702214477299*u_old[1]\
                   - 0.275607596716317*u_old[2] - 0.302855902188632*u_old[3]\
                   + 0.023993033880123*u_old[4] + 0.054021710000147*u_old[5]
        
        
        u_new[1] = - 0.460835552271137*u_old[0] - 0.029551468764759*u_old[1]\
                   + 0.393228064437077*u_old[2] + 0.105980141988699*u_old[3]\
                   + 0.030739159126096*u_old[4] - 0.039560344515975*u_old[5]
                   
                   
        u_new[2] =  0.096333139747688*u_old[0] - 0.629259919521358*u_old[1]\
                  - 0.306540940575249*u_old[2] + 1.258766829227482*u_old[3]\
                  - 0.589452745232692*u_old[4] + 0.197069524204597*u_old[5]\
                  - 0.026915887850467*u_old[6]
               
                  
        u_new[3] =  0.092596590361100*u_old[0] - 0.240558007429252*u_old[1]\
                  - 0.243684445410829*u_old[2] - 0.318292975073623*u_old[3]\
                  + 0.946962507015985*u_old[4] - 0.330443502642251*u_old[5]\
                  + 0.106765523632994*u_old[6] - 0.013345690454124*u_old[7]
                  
                  
        u_new[4] = - 0.015609194255155*u_old[0] + 0.033485250803205*u_old[1]\
                   + 0.041897570295664*u_old[2] - 0.450132219590077*u_old[3]\
                   - 0.593528379527062*u_old[4] + 1.406282074938901*u_old[5]\
                   - 0.550950133911491*u_old[6] + 0.146920035709731*u_old[7]\
                   - 0.018365004463716*u_old[8]
                   
                   
                   
        u_new[5] = - 0.016163101657400*u_old[0] + 0.040531470682718*u_old[1]\
                   - 0.005338567444210*u_old[2] - 0.025427590054054*u_old[3]\
                   - 0.357460227911152*u_old[4] - 0.572337434946968*u_old[5]\
                   + 1.313958528183954*u_old[6] - 0.492734448068983*u_old[7]\
                   + 0.131395852818395*u_old[8] - 0.016424481602299*u_old[9]
                   
                   
        u_new[m-5] =  0.032848963204599*u_old[m-7] - 0.394187558455186*u_old[m-6]\
                    - 0.572337434946968*u_old[m-5] + 1.257688453771463*u_old[m-4]\
                    - 0.406675342006739*u_old[m-3] + 0.120254802392339*u_old[m-2]\
                    - 0.054367610816546*u_old[m-1] + 0.016775726857039*u_old[m]
                    
                    
        u_new[m-4] = 0.036730008927433*u_old[m-6] - 0.399693508760138*u_old[m-5]\
                   - 0.593528379527062*u_old[m-4] + 1.303115093827634*u_old[m-3]\
                   - 0.402190050630647*u_old[m-2] + 0.047235828122731*u_old[m-1]\
                   + 0.008331008040049*u_old[m]  
                   
                   
        u_new[m-3] = - 0.020661154127888*u_old[m-5] - 0.327107204245208*u_old[m-4]\
                     - 0.318292975073623*u_old[m-3] + 0.624133691970994*u_old[m-2]\
                     + 0.118346128899150*u_old[m-1] - 0.076418487423426*u_old[m]
                     
                
        u_new[m-2] = - 0.008748664712219*u_old[m-5] + 0.061405392278187*u_old[m-4]\
                     - 0.491468255323896*u_old[m-3] - 0.306540940575249*u_old[m-2]\
                     + 0.885608401570530*u_old[m-1] - 0.140255933237354*u_old[m]
                     
                     
        u_new[m-1] = 0.029492540133095*u_old[m-5] + 0.021790841692085*u_old[m-4]\
                   - 0.215422101432619*u_old[m-3] - 0.279404147185598*u_old[m-2]\
                   - 0.029551468764759*u_old[m-1] + 0.473094335557797*u_old[m]
                   
                   
        u_new[m] = - 0.052048915548040*u_old[m-5] - 0.044953974933765*u_old[m-4]\
                   + 0.366971721882124*u_old[m-3] + 0.189297839436616*u_old[m-2]\
                   - 2.039450367044345*u_old[m-1] + 1.580183696207410*u_old[m]
                   
                   
        for i in range(6,m-5):
            
            u_new[i] = 0.033333333333333*u_old[i-2] - 0.400000000000000*u_old[i-1]\
                     - 0.583333333333333*u_old[i]   + 1.333333333333333*u_old[i+1]\
                     - 0.500000000000000*u_old[i+2] + 0.133333333333333*u_old[i+3]\
                     - 0.016666666666667*u_old[i+4]
                     
        return u_new/h
                    
                    
                  
                      
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
                        
                     
    if order == 7:
        
        u_new[0] =  -1.587120846970038 * u_old[0] + 2.068759097745568 * u_old[1] - 0.237161254615221 * u_old[2]\
                    -0.329862352927359 * u_old[3] + 0.031776313568303 * u_old[4] + 0.053609043198748 * u_old[5]
                    
        u_new[1] = - 0.463702498357543 * u_old[0] - 0.012686575154571 * u_old[1] + 0.354437950860383 * u_old[2]\
                   + 0.149830581921711 * u_old[3] + 0.006283775981431 * u_old[4] - 0.034163235251410 * u_old[5]
                   
        u_new[2] = 0.101270232427039 * u_old[0] - 0.681777727691031 * u_old[1] - 0.130372094340790 * u_old[2]\
                   + 0.980527672399361 * u_old[3] - 0.360626660648614 * u_old[4] + 0.102425925355228 * u_old[5]\
                   - 0.011447347501192 * u_old[6]
                   
        u_new[3] = 0.090637781058677 * u_old[0] - 0.220243795689184 * u_old[1] - 0.328114765554538 * u_old[2]\
                   - 0.136933409573552 * u_old[3] + 0.728101537236212 * u_old[4] - 0.181292965616422 * u_old[5]\
                   + 0.053587092315462 * u_old[6] - 0.005741474176657 * u_old[7]
                   
        u_new[4] = - 0.014132384863370 * u_old[0] + 0.024889867269617 * u_old[1] + 0.093295940554189 * u_old[2]\
                   - 0.630900759826953 * u_old[3] - 0.253703695582764 * u_old[4] + 1.044838434095626 * u_old[5]\
                   - 0.329705075321183 * u_old[6] + 0.073267794515818 * u_old[7]  - 0.007850120840981 * u_old[8]
        
        u_new[5] = - 0.016417886588667 * u_old[0] + 0.040959674614960 * u_old[1] - 0.013739135291768 * u_old[2]\
                   + 0.037096124729170 * u_old[3] - 0.551334643308827 * u_old[4] - 0.245402632435032 * u_old[5]\
                   + 0.985931565133756 * u_old[6] - 0.295779469540127 * u_old[7] + 0.065728771008917 * u_old[8]\
                   - 0.007042368322384 * u_old[9]
                   
                   
        u_new[m-5] = - 0.009389824429845 * u_old[m-8] + 0.098593156513376 * u_old[m-7] - 0.591558939080254 * u_old[m-6]\
                     - 0.245402632435032 * u_old[m-5] + 0.937327875498694 * u_old[m-4] - 0.222370039269531 * u_old[m-3]\
                     + 0.063012072625329 * u_old[m-2] - 0.046892232935029 * u_old[m-1] + 0.016680563512292 * u_old[m]
                     
        u_new[m-4] = - 0.010466827787974 * u_old[m-7] + 0.109901691773728 * u_old[m-6] - 0.614572168859249 * u_old[m-5]\
                     - 0.253703695582764 * u_old[m-4] + 0.995508274694755 * u_old[m-3] - 0.247302955053649 * u_old[m-2]\
                     + 0.009614356523787 * u_old[m-1] + 0.011021324291367 * u_old[m]
                     
                     
        u_new[m-3] = - 0.007655298902209 * u_old[m-6] + 0.030243581766320 * u_old[m-5] - 0.461432440844702 * u_old[m-4]\
                     - 0.136933409573552 * u_old[m-3] + 0.491788539658796 * u_old[m-2] + 0.167666797363523 * u_old[m-1]\
                     - 0.083677769468176 * u_old[m]
                     
        u_new[m-2] = - 0.022332921092875 * u_old[m-5] + 0.136047721252784 * u_old[m-4] - 0.654195007415723 * u_old[m-3]\
                     - 0.130372094340790 * u_old[m-2] + 0.790802931381985 * u_old[m-1] - 0.119950629785382 * u_old[m]
                     
                     
        u_new[m-1] =   0.029841082672069 * u_old[m-5] + 0.016267583768385 * u_old[m-4] - 0.196814495127562 * u_old[m-3]\
                     - 0.305572843948314 * u_old[m-2] - 0.012686575154571 * u_old[m-1] + 0.468965247789993 * u_old[m]
                     
                     
        u_new[m] =  - 0.052764835595357 * u_old[m-5] - 0.040746019354332 * u_old[m-4] + 0.357299100037565 * u_old[m-3]\
                    + 0.200227171966867 * u_old[m-2] - 2.045543388652316 * u_old[m-1] + 1.581527971597573 * u_old[m]
                     
                    
        for i in range(6, m-5):
            
            u_new[i] = - 0.009523809523810  * u_old[i-3] + 0.1 * u_old[i-2] - 0.6 * u_old[i-1] - 0.25 * u_old[i]\
                       + u_old[i+1] -0.3 * u_old[i+2] + 2/float(30) * u_old[i+3] -  0.007142857142857 * u_old[i+4]
            
            
                     
               
        
    return u_new/h
                  
                   
                   
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    

    