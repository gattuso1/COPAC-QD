do i=1,nstates-1

   do j=i,nstates-1

write(flag,'(i1,i1)') i,j
write(6,*) i, j, flag

if (flag .eq. '13') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) 

else if (flag .eq. '14') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) 

else if (flag .eq. '15') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) 

else if (flag .eq. '16') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) 

else if (flag .eq. '17') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '18') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '23') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) 

else if (flag .eq. '24') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) 

else if (flag .eq. '25') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag)

else if (flag .eq. '26') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag)

else if (flag .eq. '27') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '28') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '34') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag)

else if (flag .eq. '35') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag)

else if (flag .eq. '36') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag)

else if (flag .eq. '37') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag)

else if (flag .eq. '38') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag)

else if (flag .eq. '45') then
Ham(i,j) =D_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag)

else if (flag .eq. '46') then
Ham(i,j) =D_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aB,aA,linker,flag)

else if (flag .eq. '47') then
Ham(i,j) =D_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag)

else if (flag .eq. '48') then
Ham(i,j) =D_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag)

else if (flag .eq. '56') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag)

else if (flag .eq. '57') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) 

else if (flag .eq. '58') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '67') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '68') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag) + &
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB,linker,flag)

else if (flag .eq. '78') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) 

else if (flag .eq. '55') then
Ham(i,j) =D_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag)

else if (flag .eq. '66') then
Ham(i,j) =D_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag) + &
Dex_MCMCoff(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),&
Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker,flag)

else if (flag .eq. '77') then
Ham(i,j) =D_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) 

else if (flag .eq. '88') then
Ham(i,j) =D_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag) + &
Dex_MCMCoff(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),&
Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker,flag)

endif

   enddo

enddo

