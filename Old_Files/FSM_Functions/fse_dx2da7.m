function dx2da7 = fse_dx2da7(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1)
%FSE_DX2DA7
%    DX2DA7 = FSE_DX2DA7(F1,F2,A1,A2,A3,A4,A7,ALPHA,BETA,X1)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    23-Aug-2013 13:41:08

t1502 = 1.0./a3;
t1504 = F2.*a1;
t1503 = F1-t1504;
t1509 = alpha.*t1503;
t1505 = beta-t1509;
t1506 = exp(t1505);
t1507 = t1506+1.0;
t1508 = 1.0./t1507;
t1510 = a2-1.0;
t1511 = t1510.*t1503.*t1508;
t1512 = a7+1.0;
t1521 = t1502.*x1;
t1513 = -t1521+1.0;
t1514 = 1.0./t1512;
t1515 = t1513.^t1514;
t1516 = a3-x1;
t1535 = alpha.*t1516;
t1517 = -beta-t1535;
t1518 = exp(t1517);
t1519 = t1518+1.0;
t1520 = 1.0./t1519;
t1522 = a3.*t1512.*t1515;
t1523 = t1511+t1522;
t1536 = alpha.*t1523;
t1524 = beta-t1536;
t1525 = exp(t1524);
t1526 = t1525+1.0;
t1527 = 1.0./t1526;
t1528 = t1510.*t1502.*t1503.*t1514.*t1508;
t1529 = t1515+t1528;
t1530 = t1529.^t1512;
t1531 = a3.*t1530.*t1527;
t1532 = -a3+t1531+x1;
t1537 = t1520.*t1532.*t1508;
t1533 = t1511-t1537;
t1534 = t1502.*t1503.*x1;
t1538 = 1.0./t1512.^2;
t1539 = log(t1513);
t1540 = a3+t1537-x1;
t1545 = alpha.*t1540;
t1541 = beta-t1545;
t1542 = exp(t1541);
t1543 = t1542+1.0;
t1544 = 1.0./t1543;
t1546 = t1516.*t1508;
t1547 = t1537+t1546;
t1548 = log(t1529);
t1549 = t1530.*t1548;
t1550 = t1515.*t1538.*t1539;
t1551 = t1510.*t1502.*t1503.*t1508.*t1538;
t1552 = t1550+t1551;
t1553 = t1529.^a7;
t1566 = t1512.*t1552.*t1553;
t1554 = -t1566+t1549;
t1555 = a3.*t1527.*t1554;
t1556 = a3.*t1515;
t1567 = a3.*t1514.*t1515.*t1539;
t1557 = t1556-t1567;
t1558 = 1.0./t1526.^2;
t1559 = a3.*alpha.*t1530.*t1525.*t1557.*t1558;
t1560 = t1555+t1559;
t1561 = alpha.*t1533;
t1562 = beta+t1561;
t1563 = exp(t1562);
t1564 = t1563+1.0;
t1565 = 1.0./t1564;
t1568 = t1520.*t1560.*t1508.*t1565;
t1569 = a4.*t1502.*x1.*(1.0./2.0);
t1570 = t1534+t1569;
t1571 = t1570.*t1508;
t1572 = t1533.*t1565;
t1573 = a3.*t1508;
t1580 = t1544.*t1547;
t1574 = -t1580+t1573;
t1575 = a4.*t1502.*t1574.*(1.0./2.0);
t1576 = -t1534+t1571+t1572+t1575;
t1577 = alpha.*t1576;
t1578 = beta+t1577;
t1579 = exp(t1578);
t1581 = t1579+1.0;
t1582 = t1520.*t1560.*t1508.*t1544;
t1583 = 1.0./t1543.^2;
t1584 = alpha.*t1520.*t1542.*t1560.*t1508.*t1547.*t1583;
t1585 = t1582+t1584;
t1586 = a4.*t1502.*t1585.*(1.0./2.0);
t1587 = 1.0./t1564.^2;
t1589 = alpha.*t1520.*t1533.*t1560.*t1508.*t1563.*t1587;
t1588 = t1568+t1586-t1589;
dx2da7 = t1568-t1589-t1588./t1581+alpha.*1.0./t1581.^2.*t1576.*t1579.*t1588;
