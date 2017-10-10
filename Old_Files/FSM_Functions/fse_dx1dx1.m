function dx1dx1 = fse_dx1dx1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1)
%FSE_DX1DX1
%    DX1DX1 = FSE_DX1DX1(F1,F2,A1,A2,A3,A4,A7,ALPHA,BETA,X1)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    23-Aug-2013 13:40:42

t518 = F2.*a1;
t517 = F1-t518;
t529 = alpha.*t517;
t519 = beta-t529;
t520 = exp(t519);
t521 = t520+1.0;
t522 = 1.0./t521;
t523 = a7+1.0;
t524 = 1.0./a3;
t536 = t524.*x1;
t525 = -t536+1.0;
t526 = 1.0./t523;
t527 = t525.^t526;
t528 = a2-1.0;
t530 = a3-x1;
t549 = alpha.*t530;
t531 = -beta-t549;
t532 = exp(t531);
t533 = t532+1.0;
t534 = 1.0./t533;
t535 = t517.*t522.*t528;
t537 = a3.*t523.*t527;
t538 = t535+t537;
t545 = alpha.*t538;
t539 = beta-t545;
t540 = exp(t539);
t541 = t540+1.0;
t542 = 1.0./t541;
t543 = t517.*t522.*t524.*t526.*t528;
t544 = t527+t543;
t546 = t544.^t523;
t547 = t526-1.0;
t548 = t525.^t547;
t550 = a3.*t542.*t546;
t551 = -a3+t550+x1;
t552 = t517.*t524.*x1;
t553 = t522.*t534.*t551;
t554 = a3+t553-x1;
t560 = alpha.*t554;
t555 = beta-t560;
t556 = exp(t555);
t557 = t556+1.0;
t558 = 1.0./t557;
t559 = t517.*t524;
t561 = t544.^a7;
t562 = t542.*t548.*t561;
t563 = 1.0./t541.^2;
t564 = a3.*alpha.*t540.*t546.*t548.*t563;
t565 = t562+t564-1.0;
t566 = t522.*t534.*t565;
t567 = 1.0./t533.^2;
t568 = alpha.*t522.*t532.*t551.*t567;
t569 = t522+t566+t568;
t570 = t558.*t569;
t571 = t522.*t530;
t572 = t553+t571;
t573 = 1.0./t557.^2;
t574 = t566+t568+1.0;
t575 = alpha.*t556.*t572.*t573.*t574;
t576 = a4.*t524.*x1.*(1.0./2.0);
t577 = t552+t576;
t578 = a3.*t522;
t585 = t558.*t572;
t579 = t578-t585;
t584 = t522.*t577;
t586 = a4.*t524.*t579.*(1.0./2.0);
t580 = t552-t584-t586+x1;
t587 = alpha.*t580;
t581 = -beta-t587;
t582 = exp(t581);
t583 = t522-1.0;
t588 = t582+1.0;
t589 = a4.*t524.*(1.0./2.0);
t590 = t559+t589;
t591 = t570+t575;
t592 = t559-t522.*t590-a4.*t524.*t591.*(1.0./2.0)+1.0;
dx1dx1 = t570+t575-(t583.*t592)./t588-alpha.*t580.*t582.*t583.*1.0./t588.^2.*t592;
