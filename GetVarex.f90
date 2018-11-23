if ((flag .eq. '13') .or. (flag .eq. '14') .or. (flag .eq. '23') .or. (flag .eq. '24')) then
inA_inA='iioo' !OK
inA_inB='iiii'
inA_out='iioo'
inB_inA='oooo'
inB_inB='ooii'
inB_out='ooio'
out_inA='oooo'
out_inB='ooii'
out_out='oooo'
else if ((flag .eq. '15') .or. (flag .eq. '16') .or. (flag .eq. '25') .or. (flag .eq. '26') &
    .or. (flag .eq. '37') .or. (flag .eq. '38') .or. (flag .eq. '47') .or. (flag .eq. '48')) then
inA_inA= 'iiio'  !OK
inA_inB= 'iioi'
inA_out= 'iioo'
inB_inA= 'ooio'
inB_inB= 'oooi'
inB_out= 'oooo'
out_inA= 'ooio'
out_inB= 'oooi'
out_out= 'oooo'
else if ((flag .eq. '12') .or. (flag .eq. '13')) then
inA_inA='iiii' !OK 
inA_inB='iioo' 
inA_out='iioo' 
inB_inA='ooii' 
inB_inB='oooo'
inB_out='oooo' 
out_inA='ooii' 
out_inB='oooo' 
out_out='oooo'
else if ((flag .eq. '17') .or. (flag .eq. '18') .or. (flag .eq. '27') .or. (flag .eq. '28') &
    .or. (flag .eq. '35') .or. (flag .eq. '36') .or. (flag .eq. '45') .or. (flag .eq. '46')) then
inA_inA='iioi' !OK
inA_inB='iiio' 
inA_out='iioo' 
inB_inA='oooi' 
inB_inB='ooio' 
inB_out='oooo' 
out_inA='oooi' 
out_inB='ooio' 
out_out='oooo'
else if ((flag .eq. '55') .or. (flag .eq. '66') .or. (flag .eq. '77') .or. (flag .eq. '88') &
    .or. (flag .eq. '56') .or. (flag .eq. '78'))  then
inA_inA='ioio' !OK
inA_inB='iooi' 
inA_out='iooo' 
inB_inA='oiio' 
inB_inB='oioi' 
inB_out='oioo' 
out_inA='ooio' 
out_inB='oooi' 
out_out='oooo' 
else if ((flag .eq. '57') .or. (flag .eq. '58') .or. (flag .eq. '67') .or. (flag .eq. '68')) then
inA_inA='iooi' !OK
inA_inB='ioio' 
inA_out='iooo' 
inB_inA='oioi' 
inB_inB='oiio' 
inB_out='oioo' 
out_inA='oooi' 
out_inB='ooio' 
out_out='oooo'
endif

