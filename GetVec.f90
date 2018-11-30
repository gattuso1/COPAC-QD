if ((flag .eq. '13') .or. (flag .eq. '14') .or. (flag .eq. '23') .or. (flag .eq. '24')) then
     ra(i) = r1Anorm(i) !hole exciton 1 
     rb(j) = r2Anorm(j) !electron exciton 1
     rc(i) = r1Bnorm(i) !hole exciton 2 
     rd(j) = r2Bnorm(j) !electron exciton 2
else if ((flag .eq. '15') .or. (flag .eq. '16') .or. (flag .eq. '25') .or. (flag .eq. '26') &
    .or. (flag .eq. '37') .or. (flag .eq. '38') .or. (flag .eq. '47') .or. (flag .eq. '48')) then
     ra(i) = r1Anorm(i)
     rb(j) = r2Anorm(j)
     rc(i) = r1Anorm(i)
     rd(j) = r2Bnorm(j)
else if ((flag .eq. '12') .or. (flag .eq. '13')) then
     ra(i) = r1Anorm(i)
     rb(j) = r2Anorm(j)
     rc(i) = r1Anorm(i)
     rd(j) = r2Anorm(j)
else if ((flag .eq. '17') .or. (flag .eq. '18') .or. (flag .eq. '27') .or. (flag .eq. '28') &
    .or. (flag .eq. '35') .or. (flag .eq. '36') .or. (flag .eq. '45') .or. (flag .eq. '46')) then
     ra(i) = r1Anorm(i)
     rb(j) = r2Anorm(j)
     rc(i) = r1Bnorm(i)
     rd(j) = r2Anorm(j)
else if ((flag .eq. '55') .or. (flag .eq. '66') .or. (flag .eq. '77') .or. (flag .eq. '88') &
    .or. (flag .eq. '56') .or. (flag .eq. '78'))  then
     ra(i) = r1Anorm(i)
     rb(j) = r2Bnorm(j)
     rc(i) = r1Anorm(i)
     rd(j) = r2Bnorm(j)
else if ((flag .eq. '57') .or. (flag .eq. '58') .or. (flag .eq. '67') .or. (flag .eq. '68')) then
     ra(i) = r1Anorm(i)
     rb(j) = r2Bnorm(j)
     rc(i) = r1Bnorm(i)
     rd(j) = r2Anorm(j)
endif

