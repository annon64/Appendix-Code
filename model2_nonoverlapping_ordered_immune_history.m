%% Non-overlapping Commpartmental Formulation Incorporating Incomplete Cross-Protection
% Deterministic four-strain multilocus immune-history model

% Partial Epitope Coverage = PEC = P1 and P2, Full Epitope Coverage = FEC = P3

% Strains are defined by two loci: ax, ay, bx, by

% Immune outcomes are channeled through the following probabilities:
% P1 = locus 1 / PEC response: a or b
% P2 = locus 2 / PEC response: x or y
% P3 = locus 1 and locus 2 / FEC response: ax, ay, bx, or by

% gamma controls the degree of cross-protection between strains sharing
% one epitope. Infection risk against partially cross-protected strains is
% scaled by (1 - gamma).

% The model explicitly tracks all possible ordered immune histories,
% giving 671 epidemiological state variables.
% Setting up ODE function
function [T,Y] = model2_nonoverlapping_ordered_immune_history(beta_ax,beta_ay,beta_bx,beta_by,sigma_ax,sigma_ay,sigma_bx,sigma_by,gamma,P1,P2,P3,mu,time)

% Pack model parameters
pars = [beta_ax,beta_ay,beta_bx,beta_by,sigma_ax,sigma_ay,sigma_bx,sigma_by,gamma,P1,P2,P3,mu];

% Initial population vectors: y(1) = susceptible hosts, y(2:5) = initial infections with ax, ay, bx, by
% remaining compartments are initially empty immune-history states.
init_pop = [
  0.99799, 0.000501, 0.000502, 0.000503, 0.000504, ...
  zeros(1,666) 
];
% Specify use of ODE45 numerical solver, start at zero, steps of 3/360, ends at time
[T,Y]= ode45(@eRuations, 0:3/360:time,init_pop, ...
    odeset('nonnegative',1:671,'AbsTol', 1e-14), pars); 

end % end main function

% Define ODE eRuations
function dydt = eRuations(~,y,pars)

% Unpack Parameters back to original names 
beta_ax   = pars(1);
beta_ay   = pars(2);
beta_bx   = pars(3);
beta_by   = pars(4);
sigma_ax  = pars(5);
sigma_ay  = pars(6);
sigma_bx  = pars(7);
sigma_by  = pars(8);
gamma     = pars(9);
P1        = pars(10);
P2        = pars(11);
P3        = pars(12);
mu        = pars(13);

% Unpack epidemiological state variables.
% Variable names encode infection or immune history:
% Iax = currently infected with strain ax
% Ra  = recovered with immunity to epitope a
% Rx  = recovered with immunity to epitope x
% Rax_P3 = recovered with full strain-specific immunity to ax

% Compound names indicate accumulated immune history.

% For example:
% Iax_b = infection with ax in a host previously immune to b
% Rax_by = recovered with ax-type immunity after prior by-related history

S = y(1); Iax = y(2); Iay = y(3); Ibx = y(4); Iby = y(5);
Ra = y(6); Rb = y(7); Rx = y(8); Ry = y(9);
Rax_P3 = y(10); Ray_P3 = y(11); Rbx_P3 = y(12); Rby_P3 = y(13);

Iax_y = y(14); Iax_b = y(15); Iay_x = y(16); Iay_b = y(17);
Ibx_y = y(18); Ibx_a = y(19); Iby_x = y(20); Iby_a = y(21);
Ra_x = y(22); Ra_y = y(23); Rb_x = y(24); Rb_y = y(25);
Rx_y = y(26); Ra_b = y(27); Rx_a = y(28); Ry_a = y(29);
Rx_b = y(30); Ry_b = y(31);

Rax_y = y(32); Rax_b = y(33); Ray_x = y(34); Ray_b = y(35);
Rbx_y = y(36); Rbx_a = y(37); Rby_x = y(38); Rby_a = y(39);

Iby_a_x = y(40); Ibx_a_y = y(41); Iay_b_x = y(42); Iax_b_y = y(43);
Iby_x_a = y(44); Ibx_y_a = y(45); Iay_x_b = y(46); Iax_y_b = y(47);
Iax_ay_b = y(48); Iax_bx_y = y(49); Iay_ax_b = y(50); Iay_by_x = y(51);
Ibx_ax_y = y(52); Ibx_by_a = y(53); Iby_ay_x = y(54); Iby_bx_a = y(55);

Ra_b_y = y(56); Rx_b_y = y(57); Rax_b_y = y(58);
Ra_b_x = y(59); Ry_b_x = y(60); Ray_b_x = y(61);
Rb_a_y = y(62); Rx_a_y = y(63); Rbx_a_y = y(64);
Rb_a_x = y(65); Ry_a_x = y(66); Rby_a_x = y(67);

Ra_y_b = y(68); Rx_y_b = y(69); Rax_y_b = y(70);
Ra_x_b = y(71); Ry_x_b = y(72); Ray_x_b = y(73);
Rb_y_a = y(74); Rx_y_a = y(75); Rbx_y_a = y(76);
Rb_x_a = y(77); Ry_x_a = y(78); Rby_x_a = y(79);

Ra_ay_b = y(80); Rx_ay_b = y(81); Rax_ay_b = y(82);
Ra_bx_y = y(83); Rx_bx_y = y(84); Rax_bx_y = y(85);
Ra_ax_b = y(86); Ry_ax_b = y(87); Ray_ax_b = y(88);
Ra_by_x = y(89); Ry_by_x = y(90); Ray_by_x = y(91);
Rb_ax_y = y(92); Rx_ax_y = y(93); Rbx_ax_y = y(94);
Rb_by_a = y(95); Rx_by_a = y(96); Rbx_by_a = y(97);
Rb_ay_x = y(98); Ry_ay_x = y(99); Rby_ay_x = y(100);
Rb_bx_a = y(101); Ry_bx_a = y(102); Rby_bx_a = y(103);

Iay_ax = y(104); Ibx_ax = y(105); Iax_ay = y(106); Iby_ay = y(107);
Iax_bx = y(108); Iby_bx = y(109); Iay_by = y(110); Ibx_by = y(111);
Iax_by = y(112); Iay_bx = y(113); Ibx_ay = y(114); Iby_ax = y(115);

Ra_ax = y(116); Ry_ax = y(117); Rb_ax = y(118); Rx_ax = y(119);
Ra_ay = y(120); Rx_ay = y(121); Rb_ay = y(122); Ry_ay = y(123);
Ra_bx = y(124); Rx_bx = y(125); Rb_bx = y(126); Ry_bx = y(127);
Ra_by = y(128); Ry_by = y(129); Rb_by = y(130); Rx_by = y(131);

Rbx_ax = y(132); Ray_ax = y(133); Rbx_by = y(134); Rby_ay = y(135);
Rax_bx = y(136); Rax_ay = y(137); Rby_bx = y(138); Ray_by = y(139);
Rax_by = y(140); Ray_bx = y(141); Rby_ax = y(142); Rbx_ay = y(143);

Iax_ay_bx = y(144); Iay_ax_by = y(145); Ibx_ax_by = y(146); Iby_ay_bx = y(147);
Iay_by_ax = y(148); Ibx_by_ax = y(149); Iax_bx_ay = y(150); Iby_bx_ay = y(151);

Ra_ay_bx = y(152); Rx_ay_bx = y(153); Rax_ay_bx = y(154);
Ra_ax_by = y(155); Ry_ax_by = y(156); Ray_ax_by = y(157);
Rb_ax_by = y(158); Rx_ax_by = y(159); Rbx_ax_by = y(160);
Rb_ay_bx = y(161); Ry_ay_bx = y(162); Rby_ay_bx = y(163);
Ra_by_ax = y(164); Ry_by_ax = y(165); Ray_by_ax = y(166);
Rb_by_ax = y(167); Rx_by_ax = y(168); Rbx_by_ax = y(169);
Ra_bx_ay = y(170); Rx_bx_ay = y(171); Rax_bx_ay = y(172);
Rb_bx_ay = y(173); Ry_bx_ay = y(174); Rby_bx_ay = y(175);

Iax_b_ay_bx = y(176); Iax_y_ay_bx = y(177); Iax_by_ay_bx = y(178);
Iay_b_ax_by = y(179); Iay_x_ax_by = y(180); Iay_bx_ax_by = y(181);
Ibx_a_ax_by = y(182); Ibx_y_ax_by = y(183); Ibx_ay_ax_by = y(184);
Iby_a_ay_bx = y(185); Iby_x_ay_bx = y(186); Iby_ax_ay_bx = y(187);
Iay_b_by_ax = y(188); Iay_x_by_ax = y(189); Iay_bx_by_ax = y(190);
Ibx_a_by_ax = y(191); Ibx_y_by_ax = y(192); Ibx_ay_by_ax = y(193);
Iax_b_bx_ay = y(194); Iax_y_bx_ay = y(195); Iax_by_bx_ay = y(196);
Iby_a_bx_ay = y(197); Iby_x_bx_ay = y(198); Iby_ax_bx_ay = y(199);

Ra_b_ay_bx = y(200); Rx_b_ay_bx = y(201); Rax_b_ay_bx = y(202);
Ra_y_ay_bx = y(203); Rx_y_ay_bx = y(204); Rax_y_ay_bx = y(205);
Ra_by_ay_bx = y(206); Rx_by_ay_bx = y(207); Rax_by_ay_bx = y(208);
Ra_b_ax_by = y(209); Ry_b_ax_by = y(210); Ray_b_ax_by = y(211);
Ra_x_ax_by = y(212); Ry_x_ax_by = y(213); Ray_x_ax_by = y(214);
Ra_bx_ax_by = y(215); Ry_bx_ax_by = y(216); Ray_bx_ax_by = y(217);
Rb_a_ax_by = y(218); Rx_a_ax_by = y(219); Rbx_a_ax_by = y(220);
Rb_y_ax_by = y(221); Rx_y_ax_by = y(222); Rbx_y_ax_by = y(223);
Rb_ay_ax_by = y(224); Rx_ay_ax_by = y(225); Rbx_ay_ax_by = y(226);

Rb_a_ay_bx = y(227); Ry_a_ay_bx = y(228); Rby_a_ay_bx = y(229);
Rb_x_ay_bx = y(230); Ry_x_ay_bx = y(231); Rby_x_ay_bx = y(232);
Rb_ax_ay_bx = y(233); Ry_ax_ay_bx = y(234); Rby_ax_ay_bx = y(235);
Ra_b_bx_ay = y(236); Rx_b_bx_ay = y(237); Rax_b_bx_ay = y(238);
Ra_y_bx_ay = y(239); Rx_y_bx_ay = y(240); Rax_y_bx_ay = y(241);
Ra_by_bx_ay = y(242); Rx_by_bx_ay = y(243); Rax_by_bx_ay = y(244);
Ra_b_by_ax = y(245); Ry_b_by_ax = y(246); Ray_b_by_ax = y(247);
Ra_x_by_ax = y(248); Ry_x_by_ax = y(249); Ray_x_by_ax = y(250);
Ra_bx_by_ax = y(251); Ry_bx_by_ax = y(252); Ray_bx_by_ax = y(253);

Rb_a_by_ax = y(254); Rx_a_by_ax = y(255); Rbx_a_by_ax = y(256);
Rb_y_by_ax = y(257); Rx_y_by_ax = y(258); Rbx_y_by_ax = y(259);
Rb_ay_by_ax = y(260); Rx_ay_by_ax = y(261); Rbx_ay_by_ax = y(262);
Rb_a_bx_ay = y(263); Ry_a_bx_ay = y(264); Rby_a_bx_ay = y(265);
Rb_x_bx_ay = y(266); Ry_x_bx_ay = y(267); Rby_x_bx_ay = y(268);
Rb_ax_bx_ay = y(269); Ry_ax_bx_ay = y(270); Rby_ax_bx_ay = y(271);
Iax_b_ay = y(272); Iax_by_ay = y(273); Iax_y_bx = y(274); Iax_bx_by = y(275); 
Iay_b_ax = y(276); Iay_bx_ax = y(277); Iay_x_by = y(278); Iay_bx_by = y(279); 
Ibx_y_ax = y(280); Ibx_ay_ax = y(281); Ibx_a_by = y(282); Ibx_by_ay = y(283); 
Iay_ax_bx = y(284); Iby_ax_bx = y(285); Ibx_ax_ay = y(286); Iby_ax_ay = y(287);
Iax_by_bx = y(288); Iay_by_bx = y(289); Iax_ay_by = y(290); Ibx_ay_by = y(291);
Iby_a_bx = y(292); Iby_bx_ax = y(293); Iby_x_ay = y(294); Iby_ay_ax = y(295);
Iax_b_by = y(296); Iax_y_by = y(297); Iay_b_bx = y(298); Iay_x_bx = y(299);
Ibx_a_ay = y(300); Ibx_y_ay = y(301); Iby_a_ax = y(302); Iby_x_ax = y(303);

Ra_b_ay = y(304); Rx_b_ay = y(305); Rax_b_ay = y(306); Ra_by_ay = y(307);
Rx_by_ay = y(308); Rax_by_ay = y(309); Ra_y_bx = y(310); Rx_y_bx = y(311);
Rax_y_bx = y(312); Ra_bx_by = y(313); Rx_bx_by = y(314); Rax_bx_by = y(315);
Ra_b_ax = y(316); Ry_b_ax = y(317); Ray_b_ax = y(318); Ra_bx_ax = y(319);
Ry_bx_ax = y(320); Ray_bx_ax = y(321); Ra_x_by = y(322); Ry_x_by = y(323);
Ray_x_by = y(324); Ry_bx_by = y(325); Ray_bx_by = y(326); Rb_y_ax = y(327);
Rx_y_ax = y(328); Rbx_y_ax = y(329); Rb_ay_ax = y(330); Rx_ay_ax = y(331);
Rbx_ay_ax = y(332); Rb_a_by = y(333); Rx_a_by = y(334); Rbx_a_by = y(335);
Rb_by_ay = y(336); Rbx_by_ay = y(337); Rb_a_bx = y(338); Ry_a_bx = y(339);
Rby_a_bx = y(340); Rb_bx_ax = y(341); Rby_bx_ax = y(342); Rb_x_ay = y(343);
Ry_x_ay = y(344); Rby_x_ay = y(345); Ry_ay_ax = y(346); Rby_ay_ax = y(347);

Ra_ax_bx = y(348); Ry_ax_bx = y(349); Ray_ax_bx = y(350); Rb_ax_bx = y(351);
Rby_ax_bx = y(352); Rb_ax_ay = y(353); Rx_ax_ay = y(354); Rbx_ax_ay = y(355);
Ry_ax_ay = y(356); Rby_ax_ay = y(357); Ra_by_bx = y(358); Rx_by_bx = y(359);
Rax_by_bx = y(360); Ry_by_bx = y(361); Ray_by_bx = y(362); Ra_ay_by = y(363);
Rx_ay_by = y(364); Rax_ay_by = y(365); Rb_ay_by = y(366); Rbx_ay_by = y(367);

Ibx_a_by_ay = y(368); Ibx_ax_by_ay = y(369); Iay_x_bx_by = y(370); Iay_ax_bx_by = y(371);
Iby_a_bx_ax = y(372); Iby_ay_bx_ax = y(373); Iax_y_bx_by = y(374); Iax_ay_bx_by = y(375);
Iby_x_ay_ax = y(376); Iby_bx_ay_ax = y(377); Iax_b_by_ay = y(378); Iax_bx_by_ay = y(379);
Iay_b_bx_ax = y(380); Iay_by_bx_ax = y(381); Ibx_y_ay_ax = y(382); Ibx_by_ay_ax = y(383);
Iby_a_ax_bx = y(384); Iby_ay_ax_bx = y(385); Iay_b_ax_bx = y(386); Iay_by_ax_bx = y(387);
Iby_x_ax_ay = y(388); Iby_bx_ax_ay = y(389); Ibx_y_ax_ay = y(390); Ibx_by_ax_ay = y(391);
Iay_x_by_bx = y(392); Iay_ax_by_bx = y(393); Iax_y_by_bx = y(394); Iax_ay_by_bx = y(395);
Ibx_a_ay_by = y(396); Ibx_ax_ay_by = y(397); Iax_b_ay_by = y(398); Iax_bx_ay_by = y(399);

Rb_a_by_ay = y(400); Rx_a_by_ay = y(401); Rbx_a_by_ay = y(402); Rb_ax_by_ay = y(403);
Rx_ax_by_ay = y(404); Rbx_ax_by_ay = y(405); Ra_x_bx_by = y(406); Ry_x_bx_by = y(407);
Ray_x_bx_by = y(408); Ra_ax_bx_by = y(409); Ry_ax_bx_by = y(410); Ray_ax_bx_by = y(411);
Rb_a_bx_ax = y(412); Ry_a_bx_ax = y(413); Rby_a_bx_ax = y(414); Rb_ay_bx_ax = y(415);
Ry_ay_bx_ax = y(416); Rby_ay_bx_ax = y(417); Ra_y_bx_by = y(418); Rx_y_bx_by = y(419);
Rax_y_bx_by = y(420); Ra_ay_bx_by = y(421); Rx_ay_bx_by = y(422); Rax_ay_bx_by = y(423);
Rb_x_ay_ax = y(424); Ry_x_ay_ax = y(425); Rby_x_ay_ax = y(426); Rb_bx_ay_ax = y(427);
Ry_bx_ay_ax = y(428); Rby_bx_ay_ax = y(429);

Ra_b_by_ay = y(430); Rx_b_by_ay = y(431); Rax_b_by_ay = y(432); Ra_bx_by_ay = y(433);
Rx_bx_by_ay = y(434); Rax_bx_by_ay = y(435); Ra_b_bx_ax = y(436); Ry_b_bx_ax = y(437);
Ray_b_bx_ax = y(438); Ra_by_bx_ax = y(439); Ry_by_bx_ax = y(440); Ray_by_bx_ax = y(441);
Rb_y_ay_ax = y(442); Rx_y_ay_ax = y(443); Rbx_y_ay_ax = y(444); Rb_by_ay_ax = y(445);
Rx_by_ay_ax = y(446); Rbx_by_ay_ax = y(447); Rb_by_a_ax_bx = y(448); Ry_by_a_ax_bx = y(449);
Rby_by_a_ax_bx = y(450); Rb_by_ay_ax_bx = y(451); Ry_by_ay_ax_bx = y(452); Rby_by_ay_ax_bx = y(453);
Ra_ay_b_ax_bx = y(454); Ry_ay_b_ax_bx = y(455); Ray_ay_b_ax_bx = y(456); Ra_ay_by_ax_bx = y(457);
Ry_ay_by_ax_bx = y(458); Ray_ay_by_ax_bx = y(459);

Rb_by_x_ax_ay = y(460); Ry_by_x_ax_ay = y(461); Rby_by_x_ax_ay = y(462); Rb_by_bx_ax_ay = y(463);
Ry_by_bx_ax_ay = y(464); Rby_by_bx_ax_ay = y(465); Rb_bx_y_ax_ay = y(466); Rx_bx_y_ax_ay = y(467);
Rbx_bx_y_ax_ay = y(468); Rb_bx_by_ax_ay = y(469); Rx_bx_by_ax_ay = y(470); Rbx_bx_by_ax_ay = y(471);
Ra_ay_x_by_bx = y(472); Ry_ay_x_by_bx = y(473); Ray_ay_x_by_bx = y(474); Ra_ay_ax_by_bx = y(475);
Ry_ay_ax_by_bx = y(476); Ray_ay_ax_by_bx = y(477); Ra_ax_y_by_bx = y(478); Rx_ax_y_by_bx = y(479);
Rax_ax_y_by_bx = y(480); Ra_ax_ay_by_bx = y(481); Rx_ax_ay_by_bx = y(482); Rax_ax_ay_by_bx = y(483);
Rb_bx_a_ay_by = y(484); Rx_bx_a_ay_by = y(485); Rbx_bx_a_ay_by = y(486); Rb_bx_ax_ay_by = y(487);
Rx_bx_ax_ay_by = y(488); Rbx_bx_ax_ay_by = y(489); Ra_ax_b_ay_by = y(490); Rx_ax_b_ay_by = y(491);
Rax_ax_b_ay_by = y(492); Ra_ax_bx_ay_by = y(493); Rx_ax_bx_ay_by = y(494); Rax_ax_bx_ay_by = y(495);

Ra_b_by = y(496); Rx_b_by = y(497); Rax_b_by = y(498); Ra_y_by = y(499);
Rx_y_by = y(500); Rax_y_by = y(501); Ra_b_bx = y(502); Ry_b_bx = y(503);
Ray_b_bx = y(504); Ra_x_bx = y(505); Ry_x_bx = y(506); Ray_x_bx = y(507);
Rb_a_ay = y(508); Rx_a_ay = y(509); Rbx_a_ay = y(510); Rb_y_ay = y(511);
Rx_y_ay = y(512); Rbx_y_ay = y(513); Rb_a_ax = y(514); Ry_a_ax = y(515);
Rby_a_ax = y(516); Rb_x_ax = y(517); Ry_x_ax = y(518); Rby_x_ax = y(519);

Iay_x_b_by = y(520); Iay_ax_b_by = y(521); Ibx_a_y_by = y(522); Ibx_ax_y_by = y(523);
Iax_y_b_bx = y(524); Iax_ay_b_bx = y(525); Iby_a_x_bx = y(526); Iby_ay_x_bx = y(527);
Iby_x_a_ay = y(528); Iby_bx_a_ay = y(529); Iax_b_y_ay = y(530); Iax_bx_y_ay = y(531);
Ibx_y_a_ax = y(532); Ibx_by_a_ax = y(533); Iay_b_x_ax = y(534); Iay_by_x_ax = y(535);

Ra_x_b_by = y(536); Ry_x_b_by = y(537); Ray_x_b_by = y(538); Ra_ax_b_by = y(539);
Ry_ax_b_by = y(540); Ray_ax_b_by = y(541); Rb_a_y_by = y(542); Rx_a_y_by = y(543);
Rbx_a_y_by = y(544); Rb_ax_y_by = y(545); Rx_ax_y_by = y(546); Rbx_ax_y_by = y(547);
Ra_y_b_bx = y(548); Rx_y_b_bx = y(549); Rax_y_b_bx = y(550); Ra_ay_b_bx = y(551);
Rx_ay_b_bx = y(552); Rax_ay_b_bx = y(553); Rb_a_x_bx = y(554); Ry_a_x_bx = y(555);
Rby_a_x_bx = y(556); Rb_ay_x_bx = y(557); Ry_ay_x_bx = y(558); Rby_ay_x_bx = y(559);
Rb_x_a_ay = y(560); Ry_x_a_ay = y(561); Rby_x_a_ay = y(562); Rb_bx_a_ay = y(563);
Ry_bx_a_ay = y(564); Rby_bx_a_ay = y(565);

Ra_b_y_ay = y(566); Rx_b_y_ay = y(567); Rax_b_y_ay = y(568); Ra_bx_y_ay = y(569);
Rx_bx_y_ay = y(570); Rax_bx_y_ay = y(571); Rb_y_a_ax = y(572); Rx_y_a_ax = y(573);
Rbx_y_a_ax = y(574); Rb_by_a_ax = y(575); Rx_by_a_ax = y(576); Rbx_by_a_ax = y(577);
Ra_b_x_ax = y(578); Ry_b_x_ax = y(579); Ray_b_x_ax = y(580); Ra_by_x_ax = y(581);
Ry_by_x_ax = y(582); Ray_by_x_ax = y(583);

Iax_y_ay = y(584); Iax_b_bx = y(585); Iay_x_ax = y(586); Iay_b_by = y(587);
Ibx_a_ax = y(588); Ibx_y_by = y(589); Iby_a_ay = y(590); Iby_x_bx = y(591);
Ra_x_ax = y(592); Ra_y_ay = y(593); Rb_x_bx = y(594); Rb_y_by = y(595);
Rx_a_ax = y(596); Ry_a_ay = y(597); Rx_b_bx = y(598); Ry_b_by = y(599);
Rbx_a_ax = y(600); Ray_x_ax = y(601); Rby_a_ay = y(602); Rax_y_ay = y(603);
Rby_x_bx = y(604); Rax_b_bx = y(605); Rbx_y_by = y(606); Ray_b_by = y(607);

Iax_b_y_by = y(608); Iay_b_x_bx = y(609); Ibx_a_y_ay = y(610); Iby_a_x_ax = y(611);
Iby_x_a_ax = y(612); Ibx_y_a_ay = y(613); Iay_x_b_bx = y(614); Iax_y_b_by = y(615);
Iax_ay_b_by = y(616); Iax_bx_y_by = y(617); Iay_ax_b_bx = y(618); Iay_by_x_bx = y(619);
Ibx_ax_y_ay = y(620); Ibx_by_a_ay = y(621); Iby_ay_x_ax = y(622); Iby_bx_a_ax = y(623);

Ra_b_y_by = y(624); Rx_b_y_by = y(625); Rax_b_y_by = y(626); Ra_b_x_bx = y(627);
Ry_b_x_bx = y(628); Ray_b_x_bx = y(629); Rb_a_y_ay = y(630); Rx_a_y_ay = y(631);
Rbx_a_y_ay = y(632); Rb_a_x_ax = y(633); Ry_a_x_ax = y(634); Rby_a_x_ax = y(635);
Rb_x_a_ax = y(636); Ry_x_a_ax = y(637); Rby_x_a_ax = y(638); Rb_y_a_ay = y(639);
Rx_y_a_ay = y(640); Rbx_y_a_ay = y(641);

Ra_x_b_bx = y(642); Ry_x_b_bx = y(643); Ray_x_b_bx = y(644); Ra_y_b_by = y(645);
Rx_y_b_by = y(646); Rax_y_b_by = y(647); Ra_ay_b_by = y(648); Rx_ay_b_by = y(649);
Rax_ay_b_by = y(650); Ra_bx_y_by = y(651); Rx_bx_y_by = y(652); Rax_bx_y_by = y(653);
Ra_ax_b_bx = y(654); Ry_ax_b_bx = y(655); Ray_ax_b_bx = y(656); Ra_by_x_bx = y(657);
Ry_by_x_bx = y(658); Ray_by_x_bx = y(659); Rb_ax_y_ay = y(660); Rx_ax_y_ay = y(661);
Rbx_ax_y_ay = y(662); Rb_by_a_ay = y(663); Rx_by_a_ay = y(664); Rbx_by_a_ay = y(665);
Rb_ay_x_ax = y(666); Ry_ay_x_ax = y(667); Rby_ay_x_ax = y(668); Rb_bx_a_ax = y(669);
Ry_bx_a_ax = y(670); Rby_bx_a_ax = y(671);

% Total population size used to normalise force of infection (FOI)
N = sum(y); 

% Calculation of strain-specific FOI (lambda_ij)
% Each lambda term sums all compartments currently infectious with that strain,
% including hosts with different prior immune histories.
lambda_ax = beta_ax*( ...
    Iax + Iax_y + Iax_b + Iax_b_y + Iax_y_b + ...
    Iax_ay_b + Iax_bx_y + ...
    Iax_ay + Iax_bx + Iax_by + ...
    Iax_ay_bx + Iax_bx_ay + ...
    Iax_b_ay_bx + Iax_y_ay_bx + Iax_by_ay_bx + ...
    Iax_b_bx_ay + Iax_y_bx_ay + Iax_by_bx_ay + ...
    Iax_b_ay + Iax_by_ay + Iax_y_bx + Iax_bx_by + ...
    Iax_by_bx + Iax_ay_by + ...
    Iax_b_by + Iax_y_by + ...
    Iax_y_bx_by + Iax_ay_bx_by + ...
    Iax_b_by_ay + Iax_bx_by_ay + ...
    Iax_y_by_bx + Iax_ay_by_bx + ...
    Iax_b_ay_by + Iax_bx_ay_by + ...
    Iax_y_b_bx + Iax_ay_b_bx + ...
    Iax_b_y_ay + Iax_bx_y_ay + ...
    Iax_y_ay + Iax_b_bx + ...
    Iax_b_y_by + Iax_y_b_by + ...
    Iax_ay_b_by + Iax_bx_y_by ) / N;

lambda_ay = beta_ay*( ...
    Iay + Iay_x + Iay_b + Iay_b_x + Iay_x_b + ...
    Iay_ax_b + Iay_by_x + ...
    Iay_ax + Iay_by + Iay_bx + ...
    Iay_ax_by + Iay_by_ax + ...
    Iay_b_ax_by + Iay_x_ax_by + Iay_bx_ax_by + ...
    Iay_b_by_ax + Iay_x_by_ax + Iay_bx_by_ax + ...
    Iay_b_ax + Iay_bx_ax + Iay_x_by + Iay_bx_by + ...
    Iay_ax_bx + Iay_by_bx + ...
    Iay_b_bx + Iay_x_bx + ...
    Iay_x_bx_by + Iay_ax_bx_by + ...
    Iay_b_bx_ax + Iay_by_bx_ax + ...
    Iay_b_ax_bx + Iay_by_ax_bx + ...
    Iay_x_by_bx + Iay_ax_by_bx + ...
    Iay_x_b_by + Iay_ax_b_by + ...
    Iay_b_x_ax + Iay_by_x_ax + ...
    Iay_x_ax + Iay_b_by + ...
    Iay_b_x_bx + Iay_x_b_bx + ...
    Iay_ax_b_bx + Iay_by_x_bx ) / N;

lambda_bx = beta_bx*( ...
    Ibx + Ibx_y + Ibx_a + Ibx_a_y + Ibx_y_a + ...
    Ibx_ax_y + Ibx_by_a + ...
    Ibx_ax + Ibx_by + Ibx_ay + ...
    Ibx_ax_by + Ibx_by_ax + ...
    Ibx_a_ax_by + Ibx_y_ax_by + Ibx_ay_ax_by + ...
    Ibx_a_by_ax + Ibx_y_by_ax + Ibx_ay_by_ax + ...
    Ibx_y_ax + Ibx_ay_ax + Ibx_a_by + Ibx_by_ay + ...
    Ibx_ax_ay + Ibx_ay_by + ...
    Ibx_a_ay + Ibx_y_ay + ...
    Ibx_a_by_ay + Ibx_ax_by_ay + ...
    Ibx_y_ay_ax + Ibx_by_ay_ax + ...
    Ibx_a_ax + Ibx_y_by + ...
    Ibx_y_ax_ay + Ibx_by_ax_ay + ...
    Ibx_a_ay_by + Ibx_ax_ay_by + ...
    Ibx_a_y_by + Ibx_ax_y_by + ...
    Ibx_y_a_ax + Ibx_by_a_ax + ...
    Ibx_a_y_ay + Ibx_y_a_ay + ...
    Ibx_ax_y_ay + Ibx_by_a_ay ) / N;

lambda_by = beta_by*( ...
    Iby + Iby_x + Iby_a + Iby_a_x + Iby_x_a + ...
    Iby_ay_x + Iby_bx_a + ...
    Iby_ay + Iby_bx + Iby_ax + ...
    Iby_ay_bx + Iby_bx_ay + ...
    Iby_a_ay_bx + Iby_x_ay_bx + Iby_ax_ay_bx + ...
    Iby_a_bx_ay + Iby_x_bx_ay + Iby_ax_bx_ay + ...
    Iby_a_bx + Iby_bx_ax + Iby_x_ay + Iby_ay_ax + ...
    Iby_ax_bx + Iby_ax_ay + ...
    Iby_a_ax + Iby_x_ax + Iby_a_ay + Iby_x_bx + ...
    Iby_a_bx_ax + Iby_ay_bx_ax + ...
    Iby_x_ay_ax + Iby_bx_ay_ax + ...
    Iby_a_ax_bx + Iby_ay_ax_bx + ...
    Iby_x_ax_ay + Iby_bx_ax_ay + ...
    Iby_a_x_bx + Iby_ay_x_bx + ...
    Iby_x_a_ay + Iby_bx_a_ay + ...
    Iby_a_x_ax + Iby_x_a_ax + ...
    Iby_ay_x_ax + Iby_bx_a_ax ) / N;

%% Defining ODEs
% Susceptible dynamics:
% births enter hosts into S at rate mu*N, infections remove susceptible hosts,
% and natural mortality removes hosts at rate mu.
dS   = mu*N - S*(lambda_ax + lambda_ay + lambda_bx + lambda_by) - mu*S;

% Primary infections from the susceptible compartment.
dIax = lambda_ax*S - (sigma_ax + mu)*Iax;
dIay = lambda_ay*S - (sigma_ay + mu)*Iay;
dIbx = lambda_bx*S - (sigma_bx + mu)*Ibx;
dIby = lambda_by*S - (sigma_by + mu)*Iby;

% Recovery from primary infection.
% Hosts enter PEC/FEC immune classes according to P1, P2, and P3.
% FEC immunity is denoted R_ij_P3
dRa = sigma_ax*Iax*P1 + sigma_ay*Iay*P1 - (lambda_bx*Ra + lambda_by*Ra + mu*Ra);
dRb = sigma_bx*Ibx*P1 + sigma_by*Iby*P1 - (lambda_ax*Rb + lambda_ay*Rb + mu*Rb);
dRx = sigma_ax*Iax*P2 + sigma_bx*Ibx*P2 - (lambda_ay*Rx + lambda_by*Rx + mu*Rx); 
dRy = sigma_ay*Iay*P2 + sigma_by*Iby*P2 - (lambda_ax*Ry + lambda_bx*Ry + mu*Ry); 

dRax_P3 = sigma_ax*Iax*P3 - (lambda_by + (1 - gamma)*(lambda_ay + lambda_bx) + mu)*Rax_P3;
dRay_P3 = sigma_ay*Iay*P3 - (lambda_bx + (1 - gamma)*(lambda_ax + lambda_by) + mu)*Ray_P3;
dRbx_P3 = sigma_bx*Ibx*P3 - (lambda_ay + (1 - gamma)*(lambda_ax + lambda_by) + mu)*Rbx_P3;
dRby_P3 = sigma_by*Iby*P3 - (lambda_ax + (1 - gamma)*(lambda_ay + lambda_bx) + mu)*Rby_P3;

% Secondary infections in hosts with partial epitope-specific immunity
dIax_y = lambda_ax*Ry - (sigma_ax + mu)*Iax_y;
dIax_b = lambda_ax*Rb - (sigma_ax + mu)*Iax_b;
dIay_x = lambda_ay*Rx - (sigma_ay + mu)*Iay_x;
dIay_b = lambda_ay*Rb - (sigma_ay + mu)*Iay_b;

dIbx_y = lambda_bx*Ry - (sigma_bx + mu)*Ibx_y;
dIbx_a = lambda_bx*Ra - (sigma_bx + mu)*Ibx_a;
dIby_x = lambda_by*Rx - (sigma_by + mu)*Iby_x;
dIby_a = lambda_by*Ra - (sigma_by + mu)*Iby_a;

% Recoveries in hosts with partial epitope-specific immunity
dRa_x = sigma_ay*Iay_x*P1 - (lambda_by + mu)*Ra_x;
dRa_y = sigma_ax*Iax_y*P1 - (lambda_bx + mu)*Ra_y;
dRb_x = sigma_by*Iby_x*P1 - (lambda_ay + mu)*Rb_x;
dRb_y = sigma_bx*Ibx_y*P1 - (lambda_ax + mu)*Rb_y;
dRx_y = sigma_ax*Iax_y*P2 + sigma_ay*Iay_x*P2 + sigma_bx*Ibx_y*P2 + sigma_by*Iby_x*P2 - mu*Rx_y;
dRa_b = sigma_ax*Iax_b*P1 + sigma_ay*Iay_b*P1 + sigma_bx*Ibx_a*P1 + sigma_by*Iby_a*P1 - mu*Ra_b;

dRx_a = sigma_bx*Ibx_a*P2 - (lambda_by + mu)*Rx_a;
dRy_a = sigma_by*Iby_a*P2 - (lambda_bx + mu)*Ry_a;
dRx_b = sigma_ax*Iax_b*P2 - (lambda_ay + mu)*Rx_b;
dRy_b = sigma_ay*Iay_b*P2 - (lambda_ax + mu)*Ry_b;

dRax_y = sigma_ax*Iax_y*P3 - Rax_y*(lambda_bx + mu);
dRax_b = sigma_ax*Iax_b*P3 - Rax_b*(lambda_ay + mu);
dRay_x = sigma_ay*Iay_x*P3 - Ray_x*(lambda_by + mu);
dRay_b = sigma_ay*Iay_b*P3 - Ray_b*(lambda_ax + mu);
dRbx_y = sigma_bx*Ibx_y*P3 - Rbx_y*(lambda_ax + mu);
dRbx_a = sigma_bx*Ibx_a*P3 - Rbx_a*(lambda_by + mu);
dRby_x = sigma_by*Iby_x*P3 - Rby_x*(lambda_ay + mu);
dRby_a = sigma_by*Iby_a*P3 - Rby_a*(lambda_bx + mu);

dIby_a_x = lambda_by*Ra_x - (sigma_by + mu)*Iby_a_x;
dIbx_a_y = lambda_bx*Ra_y - (sigma_bx + mu)*Ibx_a_y;
dIay_b_x = lambda_ay*Rb_x - (sigma_ay + mu)*Iay_b_x;
dIax_b_y = lambda_ax*Rb_y - (sigma_ax + mu)*Iax_b_y;

dIby_x_a = lambda_by*Rx_a - (sigma_by + mu)*Iby_x_a;
dIbx_y_a = lambda_bx*Ry_a - (sigma_bx + mu)*Ibx_y_a;
dIay_x_b = lambda_ay*Rx_b - (sigma_ay + mu)*Iay_x_b;
dIax_y_b = lambda_ax*Ry_b - (sigma_ax + mu)*Iax_y_b;

dIax_ay_b = lambda_ax*Ray_b - Iax_ay_b*(sigma_ax + mu);
dIax_bx_y = lambda_ax*Rbx_y - Iax_bx_y*(sigma_ax + mu);
dIay_ax_b = lambda_ay*Rax_b - Iay_ax_b*(sigma_ay + mu);
dIay_by_x = lambda_ay*Rby_x - Iay_by_x*(sigma_ay + mu);
dIbx_ax_y = lambda_bx*Rax_y - Ibx_ax_y*(sigma_bx + mu);
dIbx_by_a = lambda_bx*Rby_a - Ibx_by_a*(sigma_bx + mu);
dIby_ay_x = lambda_by*Ray_x - Iby_ay_x*(sigma_by + mu);
dIby_bx_a = lambda_by*Rbx_a - Iby_bx_a*(sigma_by + mu);

dRa_b_y = sigma_ax*Iax_b_y*P1 - Ra_b_y*mu;
dRx_b_y = sigma_ax*Iax_b_y*P2 - Rx_b_y*mu;
dRax_b_y = sigma_ax*Iax_b_y*P3 - Rax_b_y*mu;

dRa_b_x = sigma_ay*Iay_b_x*P1 - Ra_b_x*mu;
dRy_b_x = sigma_ay*Iay_b_x*P2 - Ry_b_x*mu;
dRay_b_x = sigma_ay*Iay_b_x*P3 - Ray_b_x*mu;

dRb_a_y = sigma_bx*Ibx_a_y*P1 - Rb_a_y*mu;
dRx_a_y = sigma_bx*Ibx_a_y*P2 - Rx_a_y*mu;
dRbx_a_y = sigma_bx*Ibx_a_y*P3 - Rbx_a_y*mu;

dRb_a_x = sigma_by*Iby_a_x*P1 - Rb_a_x*mu;
dRy_a_x = sigma_by*Iby_a_x*P2 - Ry_a_x*mu;
dRby_a_x = sigma_by*Iby_a_x*P3 - Rby_a_x*mu;

dRa_y_b = sigma_ax*Iax_y_b*P1 - Ra_y_b*mu;
dRx_y_b = sigma_ax*Iax_y_b*P2 - Rx_y_b*mu;
dRax_y_b = sigma_ax*Iax_y_b*P3 - Rax_y_b*mu;

dRa_x_b = sigma_ay*Iay_x_b*P1 - Ra_x_b*mu;
dRy_x_b = sigma_ay*Iay_x_b*P2 - Ry_x_b*mu;
dRay_x_b = sigma_ay*Iay_x_b*P3 - Ray_x_b*mu;

dRb_y_a = sigma_bx*Ibx_y_a*P1 - Rb_y_a*mu;
dRx_y_a = sigma_bx*Ibx_y_a*P2 - Rx_y_a*mu;
dRbx_y_a = sigma_bx*Ibx_y_a*P3 - Rbx_y_a*mu;

dRb_x_a = sigma_by*Iby_x_a*P1 - Rb_x_a*mu;
dRy_x_a = sigma_by*Iby_x_a*P2 - Ry_x_a*mu;
dRby_x_a = sigma_by*Iby_x_a*P3 - Rby_x_a*mu;

dRa_ay_b   = sigma_ax*Iax_ay_b*P1 - Ra_ay_b*mu;
dRx_ay_b   = sigma_ax*Iax_ay_b*P2 - Rx_ay_b*mu;
dRax_ay_b  = sigma_ax*Iax_ay_b*P3 - Rax_ay_b*mu;

dRa_bx_y   = sigma_ax*Iax_bx_y*P1 - Ra_bx_y*mu;
dRx_bx_y   = sigma_ax*Iax_bx_y*P2 - Rx_bx_y*mu;
dRax_bx_y  = sigma_ax*Iax_bx_y*P3 - Rax_bx_y*mu;

dRa_ax_b   = sigma_ay*Iay_ax_b*P1 - Ra_ax_b*mu;
dRy_ax_b   = sigma_ay*Iay_ax_b*P2 - Ry_ax_b*mu;
dRay_ax_b  = sigma_ay*Iay_ax_b*P3 - Ray_ax_b*mu;

dRa_by_x   = sigma_ay*Iay_by_x*P1 - Ra_by_x*mu;
dRy_by_x   = sigma_ay*Iay_by_x*P2 - Ry_by_x*mu;
dRay_by_x  = sigma_ay*Iay_by_x*P3 - Ray_by_x*mu;

dRb_ax_y   = sigma_bx*Ibx_ax_y*P1 - Rb_ax_y*mu;
dRx_ax_y   = sigma_bx*Ibx_ax_y*P2 - Rx_ax_y*mu;
dRbx_ax_y  = sigma_bx*Ibx_ax_y*P3 - Rbx_ax_y*mu;

dRb_by_a   = sigma_bx*Ibx_by_a*P1 - Rb_by_a*mu;
dRx_by_a   = sigma_bx*Ibx_by_a*P2 - Rx_by_a*mu;
dRbx_by_a  = sigma_bx*Ibx_by_a*P3 - Rbx_by_a*mu;

dRb_ay_x   = sigma_by*Iby_ay_x*P1 - Rb_ay_x*mu;
dRy_ay_x   = sigma_by*Iby_ay_x*P2 - Ry_ay_x*mu;
dRby_ay_x  = sigma_by*Iby_ay_x*P3 - Rby_ay_x*mu;

dRb_bx_a   = sigma_by*Iby_bx_a*P1 - Rb_bx_a*mu;
dRy_bx_a   = sigma_by*Iby_bx_a*P2 - Ry_bx_a*mu;
dRby_bx_a  = sigma_by*Iby_bx_a*P3 - Rby_bx_a*mu;

%% Reinfection terms involving incomplete cross-protection
% Susceptibility to related strains is reduced by gamma
dIay_ax = lambda_ay*Rax_P3*(1 - gamma) - (sigma_ay + mu)*Iay_ax;
dIbx_ax = lambda_bx*Rax_P3*(1 - gamma) - (sigma_bx + mu)*Ibx_ax;

dIax_ay = lambda_ax*Ray_P3*(1 - gamma) - (sigma_ax + mu)*Iax_ay;
dIby_ay = lambda_by*Ray_P3*(1 - gamma) - (sigma_by + mu)*Iby_ay;

dIax_bx = lambda_ax*Rbx_P3*(1 - gamma) - (sigma_ax + mu)*Iax_bx;
dIby_bx = lambda_by*Rbx_P3*(1 - gamma) - (sigma_by + mu)*Iby_bx;

dIay_by = lambda_ay*Rby_P3*(1 - gamma) - (sigma_ay + mu)*Iay_by;
dIbx_by = lambda_bx*Rby_P3*(1 - gamma) - (sigma_bx + mu)*Ibx_by;

dIax_by = lambda_ax*Rby_P3 - Iax_by*(sigma_ax + mu);
dIay_bx = lambda_ay*Rbx_P3 - Iay_bx*(sigma_ay + mu);
dIbx_ay = lambda_bx*Ray_P3 - Ibx_ay*(sigma_bx + mu);
dIby_ax = lambda_by*Rax_P3 - Iby_ax*(sigma_by + mu);



dRa_ax = sigma_ay*Iay_ax*P1 - Ra_ax*(lambda_bx*(1 - gamma) + lambda_by + mu); 
dRy_ax = sigma_ay*Iay_ax*P2 + sigma_by*Iby_ax*P2 - Ry_ax*(lambda_bx*(1 - gamma) + mu);

dRb_ax = sigma_bx*Ibx_ax*P1 + sigma_by*Iby_ax*P1 - Rb_ax*(lambda_ay*(1 - gamma) + mu); 
dRx_ax = sigma_bx*Ibx_ax*P2 - Rx_ax*(lambda_ay*(1 - gamma) + lambda_by + mu); 

dRa_ay = sigma_ax*Iax_ay*P1 - Ra_ay*(lambda_by*(1 - gamma) + lambda_bx + mu); 
dRx_ay = sigma_ax*Iax_ay*P2 + sigma_bx*Ibx_ay*P2 - Rx_ay*(lambda_by*(1 - gamma) + mu); 

dRb_ay = sigma_by*Iby_ay*P1 + sigma_bx*Ibx_ay*P1 - Rb_ay*(lambda_ax*(1 - gamma) + mu); 
dRy_ay = sigma_by*Iby_ay*P2 - Ry_ay*(lambda_ax*(1 - gamma) + lambda_bx + mu); 

dRa_bx = sigma_ax*Iax_bx*P1 + sigma_ay*Iay_bx*P1 - Ra_bx*(lambda_by*(1 - gamma) + mu); 
dRx_bx = sigma_ax*Iax_bx*P2 - Rx_bx*(lambda_by*(1 - gamma) + lambda_ay + mu);

dRb_bx = sigma_by*Iby_bx*P1 - Rb_bx*(lambda_ax*(1 - gamma) + lambda_ay + mu); 
dRy_bx = sigma_by*Iby_bx*P2 + sigma_ay*Iay_bx*P2 - Ry_bx*(lambda_ax*(1 - gamma) + mu); 

dRa_by = sigma_ay*Iay_by*P1 + sigma_ax*Iax_by*P1 - Ra_by*(lambda_bx*(1 - gamma) + mu); 
dRy_by = sigma_ay*Iay_by*P2 - Ry_by*(lambda_bx*(1 - gamma) + lambda_ax + mu); 

dRb_by = sigma_bx*Ibx_by*P1 - Rb_by*(lambda_ay*(1 - gamma) + lambda_ax + mu); 
dRx_by = sigma_bx*Ibx_by*P2 + sigma_ax*Iax_by*P2 - Rx_by*(lambda_ay*(1 - gamma) + mu); 


dRbx_ax = sigma_ax*Iax_bx*P3 - Rbx_ax*(lambda_ay*(1 - gamma) + lambda_by*(1 - gamma) + mu);

dRay_ax = sigma_ay*Iay_ax*P3 - Ray_ax*(lambda_bx*(1 - gamma) + lambda_by*(1 - gamma) + mu);

dRbx_by = sigma_bx*Ibx_by*P3 - Rbx_by*(lambda_ax*(1 - gamma) + lambda_ay*(1 - gamma) + mu);

dRby_ay = sigma_by*Iby_ay*P3 - Rby_ay*(lambda_ax*(1 - gamma) + lambda_bx*(1 - gamma) + mu);

dRax_bx = sigma_bx*Ibx_ax*P3 - Rax_bx*(lambda_ay*(1 - gamma) + lambda_by*(1 - gamma) + mu);
dRax_ay = sigma_ax*Iax_ay*P3 - Rax_ay*(lambda_bx*(1 - gamma) + lambda_by*(1 - gamma) + mu);
dRby_bx = sigma_by*Iby_bx*P3 - Rby_bx*(lambda_ax*(1 - gamma) + lambda_ay*(1 - gamma) + mu);
dRay_by = sigma_ay*Iay_by*P3 - Ray_by*(lambda_ax*(1 - gamma) + lambda_bx*(1 - gamma) + mu);


dRax_by = sigma_ax*Iax_by*P3 - Rax_by*(((1 - gamma)*(lambda_ay + lambda_bx)) + mu);
dRay_bx = sigma_ay*Iay_bx*P3 - Ray_bx*(((1 - gamma)*(lambda_ax + lambda_by)) + mu);

dRby_ax = sigma_by*Iby_ax*P3 - Rby_ax*(((1 - gamma)*(lambda_ay + lambda_bx)) + mu);
dRbx_ay = sigma_bx*Ibx_ay*P3 - Rbx_ay*(((1 - gamma)*(lambda_ax + lambda_by)) + mu);


dIax_ay_bx = lambda_ax*(1 - gamma)*Ray_bx - Iax_ay_bx*(sigma_ax + mu);

dIay_ax_by = lambda_ay*(1 - gamma)*Rax_by - Iay_ax_by*(sigma_ay + mu);

dIbx_ax_by = lambda_bx*(1 - gamma)*Rax_by - Ibx_ax_by*(sigma_bx + mu);

dIby_ay_bx = lambda_by*(1 - gamma)*Ray_bx - Iby_ay_bx*(sigma_by + mu);

dIay_by_ax = lambda_ay*(1 - gamma)*Rby_ax - Iay_by_ax*(sigma_ay + mu);

dIbx_by_ax = lambda_bx*(1 - gamma)*Rby_ax - Ibx_by_ax*(sigma_bx + mu);

dIax_bx_ay = lambda_ax*(1 - gamma)*Rbx_ay - Iax_bx_ay*(sigma_ax + mu);

dIby_bx_ay = lambda_by*(1 - gamma)*Rbx_ay - Iby_bx_ay*(sigma_by + mu);


dRa_ay_bx = Iax_ay_bx*sigma_ax*P1 - Ra_ay_bx*((lambda_by*(1 - gamma)) + mu);
dRx_ay_bx = Iax_ay_bx*sigma_ax*P2 - Rx_ay_bx*((lambda_by*(1 - gamma)) + mu);
dRax_ay_bx = Iax_ay_bx*sigma_ax*P3 - Rax_ay_bx*((lambda_by*(1 - gamma)) + mu);

dRa_ax_by = Iay_ax_by*sigma_ay*P1 - Ra_ax_by*((lambda_bx*(1 - gamma)) + mu);
dRy_ax_by = Iay_ax_by*sigma_ay*P2 - Ry_ax_by*((lambda_bx*(1 - gamma)) + mu);
dRay_ax_by = Iay_ax_by*sigma_ay*P3 - Ray_ax_by*((lambda_bx*(1 - gamma)) + mu);

dRb_ax_by = Ibx_ax_by*sigma_bx*P1 - Rb_ax_by*((lambda_ay*(1 - gamma)) + mu);
dRx_ax_by = Ibx_ax_by*sigma_bx*P2 - Rx_ax_by*((lambda_ay*(1 - gamma)) + mu);
dRbx_ax_by = Ibx_ax_by*sigma_bx*P3 - Rbx_ax_by*((lambda_ay*(1 - gamma)) + mu);

dRb_ay_bx = Iby_ay_bx*sigma_by*P1 - Rb_ay_bx*((lambda_ax*(1 - gamma)) + mu);
dRy_ay_bx = Iby_ay_bx*sigma_by*P2 - Ry_ay_bx*((lambda_ax*(1 - gamma)) + mu); 
dRby_ay_bx = Iby_ay_bx*sigma_by*P3 - Rby_ay_bx*((lambda_ax*(1 - gamma)) + mu); 


dRa_by_ax = Iay_by_ax*sigma_ay*P1 - Ra_by_ax*((lambda_bx*(1 - gamma)) + mu);
dRy_by_ax = Iay_by_ax*sigma_ay*P2 - Ry_by_ax*((lambda_bx*(1 - gamma)) + mu);
dRay_by_ax = Iay_by_ax*sigma_ay*P3 - Ray_by_ax*((lambda_bx*(1 - gamma)) + mu);

dRb_by_ax = Ibx_by_ax*sigma_bx*P1 - Rb_by_ax*((lambda_ay*(1 - gamma)) + mu);
dRx_by_ax = Ibx_by_ax*sigma_bx*P2 - Rx_by_ax*((lambda_ay*(1 - gamma)) + mu);
dRbx_by_ax = Ibx_by_ax*sigma_bx*P3 - Rbx_by_ax*((lambda_ay*(1 - gamma)) + mu);


dRa_bx_ay = Iax_bx_ay*sigma_ax*P1 - Ra_bx_ay*((lambda_by*(1 - gamma)) + mu);
dRx_bx_ay = Iax_bx_ay*sigma_ax*P2 - Rx_bx_ay*((lambda_by*(1 - gamma)) + mu);
dRax_bx_ay = Iax_bx_ay*sigma_ax*P3 - Rax_bx_ay*((lambda_by*(1 - gamma)) + mu);



dRb_bx_ay = Iby_bx_ay*sigma_by*P1 - Rb_bx_ay*((lambda_ax*(1 - gamma)) + mu);
dRy_bx_ay = Iby_bx_ay*sigma_by*P2 - Ry_bx_ay*((lambda_ax*(1 - gamma)) + mu); 
dRby_bx_ay = Iby_bx_ay*sigma_by*P3 - Rby_bx_ay*((lambda_ax*(1 - gamma)) + mu); 


% 3.3 Infections from 3.2
dIax_b_ay_bx = lambda_ax*(1 - gamma)*Rb_ay_bx - Iax_b_ay_bx*(sigma_ax + mu);
dIax_y_ay_bx = lambda_ax*(1 - gamma)*Ry_ay_bx - Iax_y_ay_bx*(sigma_ax + mu);
dIax_by_ay_bx = lambda_ax*(1 - gamma)*Rby_ay_bx - Iax_by_ay_bx*(sigma_ax + mu);

dIay_b_ax_by = lambda_ay*(1 - gamma)*Rb_ax_by - Iay_b_ax_by*(sigma_ay + mu);
dIay_x_ax_by = lambda_ay*(1 - gamma)*Rx_ax_by - Iay_x_ax_by*(sigma_ay + mu);
dIay_bx_ax_by = lambda_ay*(1 - gamma)*Rbx_ax_by - Iay_bx_ax_by*(sigma_ay + mu);

dIbx_a_ax_by = lambda_bx*(1 - gamma)*Ra_ax_by - Ibx_a_ax_by*(sigma_bx + mu);
dIbx_y_ax_by = lambda_bx*(1 - gamma)*Ry_ax_by - Ibx_y_ax_by*(sigma_bx + mu);
dIbx_ay_ax_by = lambda_bx*(1 - gamma)*Ray_ax_by - Ibx_ay_ax_by*(sigma_bx + mu);

dIby_a_ay_bx = lambda_by*(1 - gamma)*Ra_ay_bx - Iby_a_ay_bx*(sigma_by + mu);
dIby_x_ay_bx = lambda_by*(1 - gamma)*Rx_ay_bx - Iby_x_ay_bx*(sigma_by + mu);
dIby_ax_ay_bx = lambda_by*(1 - gamma)*Rax_ay_bx - Iby_ax_ay_bx*(sigma_by + mu);

dIay_b_by_ax = lambda_ay*(1 - gamma)*Rb_by_ax - Iay_b_by_ax*(sigma_ay + mu);
dIay_x_by_ax = lambda_ay*(1 - gamma)*Rx_by_ax - Iay_x_by_ax*(sigma_ay + mu);
dIay_bx_by_ax = lambda_ay*(1 - gamma)*Rbx_by_ax - Iay_bx_by_ax*(sigma_ay + mu);

dIbx_a_by_ax = lambda_bx*(1 - gamma)*Ra_by_ax - Ibx_a_by_ax*(sigma_bx + mu);
dIbx_y_by_ax = lambda_bx*(1 - gamma)*Ry_by_ax - Ibx_y_by_ax*(sigma_bx + mu);
dIbx_ay_by_ax = lambda_bx*(1 - gamma)*Ray_by_ax - Ibx_ay_by_ax*(sigma_bx + mu);

dIax_b_bx_ay = lambda_ax*(1 - gamma)*Rb_bx_ay - Iax_b_bx_ay*(sigma_ax + mu);
dIax_y_bx_ay = lambda_ax*(1 - gamma)*Ry_bx_ay - Iax_y_bx_ay*(sigma_ax + mu);
dIax_by_bx_ay = lambda_ax*(1 - gamma)*Rby_bx_ay - Iax_by_bx_ay*(sigma_ax + mu);

dIby_a_bx_ay = lambda_by*(1 - gamma)*Ra_bx_ay - Iby_a_bx_ay*(sigma_by + mu);
dIby_x_bx_ay = lambda_by*(1 - gamma)*Rx_bx_ay - Iby_x_bx_ay*(sigma_by + mu);
dIby_ax_bx_ay = lambda_by*(1 - gamma)*Rax_bx_ay - Iby_ax_bx_ay*(sigma_by + mu);

dRa_b_ay_bx  = sigma_ax*Iax_b_ay_bx*P1 - Ra_b_ay_bx*mu;
dRx_b_ay_bx  = sigma_ax*Iax_b_ay_bx*P2 - Rx_b_ay_bx*mu;
dRax_b_ay_bx = sigma_ax*Iax_b_ay_bx*P3 - Rax_b_ay_bx*mu;

dRa_y_ay_bx  = sigma_ax*Iax_y_ay_bx*P1 - Ra_y_ay_bx*mu;
dRx_y_ay_bx  = sigma_ax*Iax_y_ay_bx*P2 - Rx_y_ay_bx*mu;
dRax_y_ay_bx = sigma_ax*Iax_y_ay_bx*P3 - Rax_y_ay_bx*mu;

dRa_by_ay_bx  = sigma_ax*Iax_by_ay_bx*P1 - Ra_by_ay_bx*mu;
dRx_by_ay_bx  = sigma_ax*Iax_by_ay_bx*P2 - Rx_by_ay_bx*mu;
dRax_by_ay_bx = sigma_ax*Iax_by_ay_bx*P3 - Rax_by_ay_bx*mu;

dRa_b_ax_by  = sigma_ay*Iay_b_ax_by*P1 - Ra_b_ax_by*mu;
dRy_b_ax_by  = sigma_ay*Iay_b_ax_by*P2 - Ry_b_ax_by*mu;
dRay_b_ax_by = sigma_ay*Iay_b_ax_by*P3 - Ray_b_ax_by*mu;

dRa_x_ax_by  = sigma_ay*Iay_x_ax_by*P1 - Ra_x_ax_by*mu;
dRy_x_ax_by  = sigma_ay*Iay_x_ax_by*P2 - Ry_x_ax_by*mu;
dRay_x_ax_by = sigma_ay*Iay_x_ax_by*P3 - Ray_x_ax_by*mu;

dRa_bx_ax_by  = sigma_ay*Iay_bx_ax_by*P1 - Ra_bx_ax_by*mu;
dRy_bx_ax_by  = sigma_ay*Iay_bx_ax_by*P2 - Ry_bx_ax_by*mu;
dRay_bx_ax_by = sigma_ay*Iay_bx_ax_by*P3 - Ray_bx_ax_by*mu;

dRb_a_ax_by  = sigma_bx*Ibx_a_ax_by*P1 - Rb_a_ax_by*mu;
dRx_a_ax_by  = sigma_bx*Ibx_a_ax_by*P2 - Rx_a_ax_by*mu;
dRbx_a_ax_by = sigma_bx*Ibx_a_ax_by*P3 - Rbx_a_ax_by*mu;

dRb_y_ax_by  = sigma_bx*Ibx_y_ax_by*P1 - Rb_y_ax_by*mu;
dRx_y_ax_by  = sigma_bx*Ibx_y_ax_by*P2 - Rx_y_ax_by*mu;
dRbx_y_ax_by = sigma_bx*Ibx_y_ax_by*P3 - Rbx_y_ax_by*mu;

dRb_ay_ax_by  = sigma_bx*Ibx_ay_ax_by*P1 - Rb_ay_ax_by*mu;
dRx_ay_ax_by  = sigma_bx*Ibx_ay_ax_by*P2 - Rx_ay_ax_by*mu;
dRbx_ay_ax_by = sigma_bx*Ibx_ay_ax_by*P3 - Rbx_ay_ax_by*mu;

dRb_a_ay_bx  = sigma_by*Iby_a_ay_bx*P1 - Rb_a_ay_bx*mu;
dRy_a_ay_bx  = sigma_by*Iby_a_ay_bx*P2 - Ry_a_ay_bx*mu;
dRby_a_ay_bx = sigma_by*Iby_a_ay_bx*P3 - Rby_a_ay_bx*mu;

dRb_x_ay_bx  = sigma_by*Iby_x_ay_bx*P1 - Rb_x_ay_bx*mu;
dRy_x_ay_bx  = sigma_by*Iby_x_ay_bx*P2 - Ry_x_ay_bx*mu;
dRby_x_ay_bx = sigma_by*Iby_x_ay_bx*P3 - Rby_x_ay_bx*mu;

dRb_ax_ay_bx  = sigma_by*Iby_ax_ay_bx*P1 - Rb_ax_ay_bx*mu;
dRy_ax_ay_bx  = sigma_by*Iby_ax_ay_bx*P2 - Ry_ax_ay_bx*mu;
dRby_ax_ay_bx = sigma_by*Iby_ax_ay_bx*P3 - Rby_ax_ay_bx*mu;

dRa_b_bx_ay  = sigma_ax*Iax_b_bx_ay*P1 - Ra_b_bx_ay*mu;
dRx_b_bx_ay  = sigma_ax*Iax_b_bx_ay*P2 - Rx_b_bx_ay*mu;
dRax_b_bx_ay = sigma_ax*Iax_b_bx_ay*P3 - Rax_b_bx_ay*mu;

dRa_y_bx_ay  = sigma_ax*Iax_y_bx_ay*P1 - Ra_y_bx_ay*mu;
dRx_y_bx_ay  = sigma_ax*Iax_y_bx_ay*P2 - Rx_y_bx_ay*mu;
dRax_y_bx_ay = sigma_ax*Iax_y_bx_ay*P3 - Rax_y_bx_ay*mu;

dRa_by_bx_ay  = sigma_ax*Iax_by_bx_ay*P1 - Ra_by_bx_ay*mu;
dRx_by_bx_ay  = sigma_ax*Iax_by_bx_ay*P2 - Rx_by_bx_ay*mu;
dRax_by_bx_ay = sigma_ax*Iax_by_bx_ay*P3 - Rax_by_bx_ay*mu;

dRa_b_by_ax  = sigma_ay*Iay_b_by_ax*P1 - Ra_b_by_ax*mu;
dRy_b_by_ax  = sigma_ay*Iay_b_by_ax*P2 - Ry_b_by_ax*mu;
dRay_b_by_ax = sigma_ay*Iay_b_by_ax*P3 - Ray_b_by_ax*mu;

dRa_x_by_ax  = sigma_ay*Iay_x_by_ax*P1 - Ra_x_by_ax*mu;
dRy_x_by_ax  = sigma_ay*Iay_x_by_ax*P2 - Ry_x_by_ax*mu;
dRay_x_by_ax = sigma_ay*Iay_x_by_ax*P3 - Ray_x_by_ax*mu;

dRa_bx_by_ax  = sigma_ay*Iay_bx_by_ax*P1 - Ra_bx_by_ax*mu;
dRy_bx_by_ax  = sigma_ay*Iay_bx_by_ax*P2 - Ry_bx_by_ax*mu;
dRay_bx_by_ax = sigma_ay*Iay_bx_by_ax*P3 - Ray_bx_by_ax*mu;


dRb_a_by_ax  = sigma_bx*Ibx_a_by_ax*P1 - Rb_a_by_ax*mu;
dRx_a_by_ax  = sigma_bx*Ibx_a_by_ax*P2 - Rx_a_by_ax*mu;
dRbx_a_by_ax = sigma_bx*Ibx_a_by_ax*P3 - Rbx_a_by_ax*mu;

dRb_y_by_ax  = sigma_bx*Ibx_y_by_ax*P1 - Rb_y_by_ax*mu;
dRx_y_by_ax  = sigma_bx*Ibx_y_by_ax*P2 - Rx_y_by_ax*mu;
dRbx_y_by_ax = sigma_bx*Ibx_y_by_ax*P3 - Rbx_y_by_ax*mu;

dRb_ay_by_ax  = sigma_bx*Ibx_ay_by_ax*P1 - Rb_ay_by_ax*mu;
dRx_ay_by_ax  = sigma_bx*Ibx_ay_by_ax*P2 - Rx_ay_by_ax*mu;
dRbx_ay_by_ax = sigma_bx*Ibx_ay_by_ax*P3 - Rbx_ay_by_ax*mu;

dRb_a_bx_ay  = sigma_by*Iby_a_bx_ay*P1 - Rb_a_bx_ay*mu;
dRy_a_bx_ay  = sigma_by*Iby_a_bx_ay*P2 - Ry_a_bx_ay*mu;
dRby_a_bx_ay = sigma_by*Iby_a_bx_ay*P3 - Rby_a_bx_ay*mu;

dRb_x_bx_ay  = sigma_by*Iby_x_bx_ay*P1 - Rb_x_bx_ay*mu;
dRy_x_bx_ay  = sigma_by*Iby_x_bx_ay*P2 - Ry_x_bx_ay*mu;
dRby_x_bx_ay = sigma_by*Iby_x_bx_ay*P3 - Rby_x_bx_ay*mu;

dRb_ax_bx_ay  = sigma_by*Iby_ax_bx_ay*P1 - Rb_ax_bx_ay*mu;
dRy_ax_bx_ay  = sigma_by*Iby_ax_bx_ay*P2 - Ry_ax_bx_ay*mu;
dRby_ax_bx_ay = sigma_by*Iby_ax_bx_ay*P3 - Rby_ax_bx_ay*mu;

dIax_b_ay = lambda_ax*(1 - gamma)*Rb_ay - Iax_b_ay*(sigma_ax + mu);
dIax_by_ay = lambda_ax*(1 - gamma)*Rby_ay - Iax_by_ay*(sigma_ax + mu);
dIax_y_bx = lambda_ax*(1 - gamma)* Ry_bx - Iax_y_bx*(sigma_ax + mu);
dIax_bx_by = lambda_ax*(1 - gamma)* Rbx_by - Iax_bx_by*(sigma_ax + mu);

dIay_b_ax = lambda_ay*(1 - gamma)*Rb_ax- Iay_b_ax*(sigma_ay + mu);
dIay_bx_ax = lambda_ay*(1 - gamma)*Rbx_ax - Iay_bx_ax*(sigma_ay + mu);
dIay_x_by = lambda_ay*(1 - gamma)*Rx_by - Iay_x_by*(sigma_ay + mu);
dIay_bx_by = lambda_ay*(1 - gamma)*Rbx_by - Iay_bx_by*(sigma_ay + mu);

dIbx_y_ax = lambda_bx*(1 - gamma)*Ry_ax - Ibx_y_ax*(sigma_bx + mu);
dIbx_ay_ax = lambda_bx*(1 - gamma)*Ray_ax - Ibx_ay_ax*(sigma_bx + mu);
dIbx_a_by = lambda_bx*(1 - gamma)*Ra_by - Ibx_a_by*(sigma_bx + mu);
dIbx_by_ay = lambda_bx*(1 - gamma)*Rby_ay - Ibx_by_ay*(sigma_bx + mu);

dIay_ax_bx = lambda_ay*(1 - gamma)*Rax_bx - Iay_ax_bx*(sigma_ay + mu);
dIby_ax_bx = lambda_by*(1 - gamma)*Rax_bx - Iby_ax_bx*(sigma_by + mu);
dIbx_ax_ay = lambda_bx*(1 - gamma)*Rax_ay - Ibx_ax_ay*(sigma_bx + mu);
dIby_ax_ay = lambda_by*(1 - gamma)*Rax_ay - Iby_ax_ay*(sigma_by + mu);
dIax_by_bx = lambda_ax*(1 - gamma)*Rby_bx - Iax_by_bx*(sigma_ax + mu);
dIay_by_bx = lambda_ay*(1 - gamma)*Rby_bx - Iay_by_bx*(sigma_ay + mu);

dIax_ay_by = lambda_ax*(1 - gamma)*Ray_by - Iax_ay_by*(sigma_ax + mu);
dIbx_ay_by = lambda_bx*(1 - gamma)*Ray_by - Ibx_ay_by*(sigma_bx + mu);

dIby_a_bx = lambda_by*(1 - gamma)*Ra_bx - Iby_a_bx*(sigma_by + mu);
dIby_bx_ax = lambda_by*(1 - gamma)*Rbx_ax - Iby_bx_ax*(sigma_by + mu);
dIby_x_ay = lambda_by*(1 - gamma)*Rx_ay - Iby_x_ay*(sigma_by + mu);
dIby_ay_ax = lambda_by*(1 - gamma)*Ray_ax - Iby_ay_ax*(sigma_by + mu);

dIax_b_by = lambda_ax*Rb_by - Iax_b_by*(sigma_ax + mu);
dIax_y_by = lambda_ax*Ry_by - Iax_y_by*(sigma_ax + mu);

dIay_b_bx = lambda_ay*Rb_bx - Iay_b_bx*(sigma_ay + mu);
dIay_x_bx = lambda_ay*Rx_bx - Iay_x_bx*(sigma_ay + mu);

dIbx_a_ay = lambda_bx*Ra_ay - Ibx_a_ay*(sigma_bx + mu);
dIbx_y_ay = lambda_bx*Ry_ay - Ibx_y_ay*(sigma_bx + mu);

dIby_a_ax = lambda_by*Ra_ax - Iby_a_ax*(sigma_by + mu);
dIby_x_ax = lambda_by*Rx_ax - Iby_x_ax*(sigma_by + mu);

dRa_b_ay = sigma_ax*Iax_b_ay*P1 - Ra_b_ay*mu;
dRx_b_ay = sigma_ax*Iax_b_ay*P2 - Rx_b_ay*mu;
dRax_b_ay = sigma_ax*Iax_b_ay*P3 - Rax_b_ay*mu;

dRa_by_ay = sigma_ax*Iax_by_ay*P1 - Ra_by_ay*(lambda_bx*(1 - gamma) + mu);
dRx_by_ay = sigma_ax*Iax_by_ay*P2 + sigma_bx*Ibx_by_ay*P2 - Rx_by_ay*mu;
dRax_by_ay = sigma_ax*Iax_by_ay*P3 - Rax_by_ay*(lambda_bx*(1 - gamma) + mu); 

dRa_y_bx = sigma_ax*Iax_y_bx*P1 - Ra_y_bx*mu;
dRx_y_bx = sigma_ax*Iax_y_bx*P2 - Rx_y_bx*mu;
dRax_y_bx = sigma_ax*Iax_y_bx*P3 - Rax_y_bx*mu;

dRa_bx_by = sigma_ax*Iax_bx_by*P1 + sigma_ay*Iay_bx_by*P1 - Ra_bx_by*mu; 
dRx_bx_by = sigma_ax*Iax_bx_by*P2 - Rx_bx_by*(lambda_ay*(1 - gamma) + mu);
dRax_bx_by = sigma_ax*Iax_bx_by*P3 - Rax_bx_by*(lambda_ay*(1 - gamma) + mu);

dRa_b_ax = sigma_ay*Iay_b_ax*P1 - Ra_b_ax*mu;
dRy_b_ax = sigma_ay*Iay_b_ax*P2 - Ry_b_ax*mu;
dRay_b_ax = sigma_ay*Iay_b_ax*P3 - Ray_b_ax*mu;

dRa_bx_ax = sigma_ay*Iay_bx_ax*P1 - Ra_bx_ax*(lambda_by*(1 - gamma) + mu);
dRy_bx_ax = sigma_ay*Iay_bx_ax*P2 + sigma_by*Iby_bx_ax*P2 - Ry_bx_ax*mu;
dRay_bx_ax = sigma_ay*Iay_bx_ax*P3 - Ray_bx_ax*(lambda_by*(1 - gamma) + mu);

dRa_x_by = sigma_ay*Iay_x_by*P1 - Ra_x_by*mu;
dRy_x_by = sigma_ay*Iay_x_by*P2 - Ry_x_by*mu;
dRay_x_by = sigma_ay*Iay_x_by*P3 - Ray_x_by*mu;

dRy_bx_by = sigma_ay*Iay_bx_by*P2 - Ry_bx_by*(lambda_ax*(1 - gamma) + mu);
dRay_bx_by = sigma_ay*Iay_bx_by*P3 - Ray_bx_by*(lambda_ax*(1 - gamma) + mu);

dRb_y_ax = sigma_bx*Ibx_y_ax*P1 - Rb_y_ax*mu; 
dRx_y_ax = sigma_bx*Ibx_y_ax*P2 - Rx_y_ax*mu;
dRbx_y_ax = sigma_bx*Ibx_y_ax*P3 - Rbx_y_ax*mu;

dRb_ay_ax = sigma_bx*Ibx_ay_ax*P1 + sigma_by*Iby_ay_ax*P1 - Rb_ay_ax*mu;
dRx_ay_ax = sigma_bx*Ibx_ay_ax*P2 - Rx_ay_ax*(lambda_by*(1 - gamma) + mu);
dRbx_ay_ax = sigma_bx*Ibx_ay_ax*P3 - Rbx_ay_ax*(lambda_by*(1 - gamma) + mu);

dRb_a_by = sigma_bx*Ibx_a_by*P1 - Rb_a_by*mu;
dRx_a_by = sigma_bx*Ibx_a_by*P2 - Rx_a_by*mu;
dRbx_a_by = sigma_bx*Ibx_a_by*P3 - Rbx_a_by*mu;

dRb_by_ay = sigma_bx*Ibx_by_ay*P1 - Rb_by_ay*(lambda_ax*(1 - gamma) + mu);
dRbx_by_ay = sigma_bx*Ibx_by_ay*P3 - Rbx_by_ay*(lambda_ax*(1 - gamma) + mu);

dRb_a_bx = sigma_by*Iby_a_bx*P1 - Rb_a_bx*mu;
dRy_a_bx = sigma_by*Iby_a_bx*P2 - Ry_a_bx*mu;
dRby_a_bx = sigma_by*Iby_a_bx*P3 - Rby_a_bx*mu;

dRb_bx_ax = sigma_by*Iby_bx_ax*P1 - Rb_bx_ax*(lambda_ay*(1 - gamma) + mu);
dRby_bx_ax = sigma_by*Iby_bx_ax*P3 - Rby_bx_ax*(lambda_ay*(1 - gamma) + mu);

dRb_x_ay = sigma_by*Iby_x_ay*P1 - Rb_x_ay*mu;
dRy_x_ay = sigma_by*Iby_x_ay*P2 - Ry_x_ay*mu;
dRby_x_ay = sigma_by*Iby_x_ay*P3 - Rby_x_ay*mu;

dRy_ay_ax = sigma_by*Iby_ay_ax*P2 - Ry_ay_ax*(lambda_bx*(1 - gamma) + mu);
dRby_ay_ax = sigma_by*Iby_ay_ax*P3 - Rby_ay_ax*(lambda_bx*(1 - gamma) + mu);

dRa_ax_bx = sigma_ay*Iay_ax_bx*P1 - Ra_ax_bx*(lambda_by*(1 - gamma) + mu);
dRy_ax_bx = sigma_ay*Iay_ax_bx*P2 + sigma_by*Iby_ax_bx*P2 - Ry_ax_bx*mu;
dRay_ax_bx = sigma_ay*Iay_ax_bx*P3 - Ray_ax_bx*(lambda_by*(1 - gamma) + mu);
dRb_ax_bx = sigma_by*Iby_ax_bx*P1 - Rb_ax_bx*(lambda_ay*(1 - gamma) + mu);
dRby_ax_bx = sigma_by*Iby_ax_bx*P3 - Rby_ax_bx*(lambda_ay*(1 - gamma) + mu);

dRb_ax_ay = sigma_bx*Ibx_ax_ay*P1 + sigma_by*Iby_ax_ay*P1 - Rb_ax_ay*mu;
dRx_ax_ay = sigma_bx*Ibx_ax_ay*P2 - Rx_ax_ay*(lambda_by*(1 - gamma) + mu);
dRbx_ax_ay = sigma_bx*Ibx_ax_ay*P3 - Rbx_ax_ay*(lambda_by*(1 - gamma) + mu);
dRy_ax_ay = sigma_by*Iby_ax_ay*P2 - Ry_ax_ay*(lambda_bx*(1 - gamma) + mu);
dRby_ax_ay = sigma_by*Iby_ax_ay*P3 - Rby_ax_ay*(lambda_bx*(1 - gamma) + mu);

dRa_by_bx = sigma_ax*Iax_by_bx*P1 + sigma_ay*Iay_by_bx*P1 - Ra_by_bx*mu; 
dRx_by_bx = sigma_ax*Iax_by_bx*P2 - Rx_by_bx*(lambda_ay*(1 - gamma) + mu);
dRax_by_bx = sigma_ax*Iax_by_bx*P3 - Rax_by_bx*(lambda_ay*(1 - gamma) + mu);

dRy_by_bx = sigma_ay*Iay_by_bx*P2 - Ry_by_bx*(lambda_ax*(1 - gamma) + mu);
dRay_by_bx = sigma_ay*Iay_by_bx*P3 - Ray_by_bx*(lambda_ax*(1 - gamma) + mu);

dRa_ay_by = sigma_ax*Iax_ay_by*P1 - Ra_ay_by*(lambda_bx*(1 - gamma) + mu);
dRx_ay_by = sigma_ax*Iax_ay_by*P2 + sigma_bx*Ibx_ay_by*P2 - Rx_ay_by*mu;
dRax_ay_by = sigma_ax*Iax_ay_by*P3 - Rax_ay_by*(lambda_bx*(1 - gamma) + mu); 
dRb_ay_by = sigma_bx*Ibx_ay_by*P1 - Rb_ay_by*(lambda_ax*(1 - gamma) + mu);
dRbx_ay_by = sigma_bx*Ibx_ay_by*P3 - Rbx_ay_by*(lambda_ax*(1 - gamma) + mu);

dIbx_a_by_ay = lambda_bx*(1 - gamma)*Ra_by_ay - Ibx_a_by_ay*(sigma_bx + mu);
dIbx_ax_by_ay = lambda_bx*(1 - gamma)*Rax_by_ay - Ibx_ax_by_ay*(sigma_bx + mu);

dIay_x_bx_by = lambda_ay*(1 - gamma)*Rx_bx_by - Iay_x_bx_by*(sigma_ay + mu);
dIay_ax_bx_by = lambda_ay*(1 - gamma)*Rax_bx_by - Iay_ax_bx_by*(sigma_ay + mu);

dIby_a_bx_ax = lambda_by*(1 - gamma)*Ra_bx_ax - Iby_a_bx_ax*(sigma_by + mu);
dIby_ay_bx_ax = lambda_by*(1 - gamma)*Ray_bx_ax - Iby_ay_bx_ax*(sigma_by + mu);

dIax_y_bx_by = lambda_ax*(1 - gamma)*Ry_bx_by - Iax_y_bx_by*(sigma_ax + mu);
dIax_ay_bx_by = lambda_ax*(1 - gamma)*Ray_bx_by - Iax_ay_bx_by*(sigma_ax + mu);

dIby_x_ay_ax = lambda_by*(1 - gamma)*Rx_ay_ax - Iby_x_ay_ax*(sigma_by + mu);
dIby_bx_ay_ax = lambda_by*(1 - gamma)*Rbx_ay_ax - Iby_bx_ay_ax*(sigma_by + mu);

dIax_b_by_ay = lambda_ax*(1 - gamma)*Rb_by_ay - Iax_b_by_ay*(sigma_ax + mu);
dIax_bx_by_ay = lambda_ax*(1 - gamma)*Rbx_by_ay - Iax_bx_by_ay*(sigma_ax + mu);

dIay_b_bx_ax = lambda_ay*(1 - gamma)*Rb_bx_ax - Iay_b_bx_ax*(sigma_ay + mu);
dIay_by_bx_ax = lambda_ay*(1 - gamma)*Rby_bx_ax - Iay_by_bx_ax*(sigma_ay + mu);

dIbx_y_ay_ax = lambda_bx*(1 - gamma)*Ry_ay_ax - Ibx_y_ay_ax*(sigma_bx + mu);
dIbx_by_ay_ax = lambda_bx*(1 - gamma)*Rby_ay_ax - Ibx_by_ay_ax*(sigma_bx + mu);



dIby_a_ax_bx  = lambda_by*(1 - gamma)*Ra_ax_bx  - (sigma_by + mu)*Iby_a_ax_bx;
dIby_ay_ax_bx = lambda_by*(1 - gamma)*Ray_ax_bx - (sigma_by + mu)*Iby_ay_ax_bx;

dIay_b_ax_bx  = lambda_ay*(1 - gamma)*Rb_ax_bx  - (sigma_ay + mu)*Iay_b_ax_bx;
dIay_by_ax_bx = lambda_ay*(1 - gamma)*Rby_ax_bx - (sigma_ay + mu)*Iay_by_ax_bx;

dIby_x_ax_ay  = lambda_by*(1 - gamma)*Rx_ax_ay  - (sigma_by + mu)*Iby_x_ax_ay;
dIby_bx_ax_ay = lambda_by*(1 - gamma)*Rbx_ax_ay - (sigma_by + mu)*Iby_bx_ax_ay;

dIbx_y_ax_ay  = lambda_bx*(1 - gamma)*Ry_ax_ay  - (sigma_bx + mu)*Ibx_y_ax_ay;
dIbx_by_ax_ay = lambda_bx*(1 - gamma)*Rby_ax_ay - (sigma_bx + mu)*Ibx_by_ax_ay;


dIay_x_by_bx  = lambda_ay*(1 - gamma)*Rx_by_bx  - (sigma_ay + mu)*Iay_x_by_bx;
dIay_ax_by_bx = lambda_ay*(1 - gamma)*Rax_by_bx - (sigma_ay + mu)*Iay_ax_by_bx;

dIax_y_by_bx  = lambda_ax*(1 - gamma)*Ry_by_bx  - (sigma_ax + mu)*Iax_y_by_bx;
dIax_ay_by_bx = lambda_ax*(1 - gamma)*Ray_by_bx - (sigma_ax + mu)*Iax_ay_by_bx;


dIbx_a_ay_by  = lambda_bx*(1 - gamma)*Ra_ay_by  - (sigma_bx + mu)*Ibx_a_ay_by;
dIbx_ax_ay_by = lambda_bx*(1 - gamma)*Rax_ay_by - (sigma_bx + mu)*Ibx_ax_ay_by;

dIax_b_ay_by  = lambda_ax*(1 - gamma)*Rb_ay_by  - (sigma_ax + mu)*Iax_b_ay_by;
dIax_bx_ay_by = lambda_ax*(1 - gamma)*Rbx_ay_by - (sigma_ax + mu)*Iax_bx_ay_by;


dRb_a_by_ay = sigma_bx*Ibx_a_by_ay*P1 - Rb_a_by_ay*mu;
dRx_a_by_ay = sigma_bx*Ibx_a_by_ay*P2 - Rx_a_by_ay*mu;
dRbx_a_by_ay = sigma_bx*Ibx_a_by_ay*P3 - Rbx_a_by_ay*mu;

dRb_ax_by_ay = sigma_bx*Ibx_ax_by_ay*P1 - Rb_ax_by_ay*mu;
dRx_ax_by_ay = sigma_bx*Ibx_ax_by_ay*P2 - Rx_ax_by_ay*mu;
dRbx_ax_by_ay = sigma_bx*Ibx_ax_by_ay*P3 - Rbx_ax_by_ay*mu;

dRa_x_bx_by = sigma_ay*Iay_x_bx_by*P1 - Ra_x_bx_by*mu;
dRy_x_bx_by = sigma_ay*Iay_x_bx_by*P2 - Ry_x_bx_by*mu;
dRay_x_bx_by = sigma_ay*Iay_x_bx_by*P3 - Ray_x_bx_by*mu;

dRa_ax_bx_by = sigma_ay*Iay_ax_bx_by*P1 - Ra_ax_bx_by*mu;
dRy_ax_bx_by = sigma_ay*Iay_ax_bx_by*P2 - Ry_ax_bx_by*mu;
dRay_ax_bx_by = sigma_ay*Iay_ax_bx_by*P3 - Ray_ax_bx_by*mu;

dRb_a_bx_ax = sigma_by*Iby_a_bx_ax*P1 - Rb_a_bx_ax*mu;
dRy_a_bx_ax = sigma_by*Iby_a_bx_ax*P2 - Ry_a_bx_ax*mu;
dRby_a_bx_ax = sigma_by*Iby_a_bx_ax*P3 - Rby_a_bx_ax*mu;

dRb_ay_bx_ax = sigma_by*Iby_ay_bx_ax*P1 - Rb_ay_bx_ax*mu;
dRy_ay_bx_ax = sigma_by*Iby_ay_bx_ax*P2 - Ry_ay_bx_ax*mu;
dRby_ay_bx_ax = sigma_by*Iby_ay_bx_ax*P3 - Rby_ay_bx_ax*mu;

dRa_y_bx_by = sigma_ax*Iax_y_bx_by*P1 - Ra_y_bx_by*mu;
dRx_y_bx_by = sigma_ax*Iax_y_bx_by*P2 - Rx_y_bx_by*mu;
dRax_y_bx_by = sigma_ax*Iax_y_bx_by*P3 - Rax_y_bx_by*mu;

dRa_ay_bx_by = sigma_ax*Iax_ay_bx_by*P1 - Ra_ay_bx_by*mu;
dRx_ay_bx_by = sigma_ax*Iax_ay_bx_by*P2 - Rx_ay_bx_by*mu;
dRax_ay_bx_by = sigma_ax*Iax_ay_bx_by*P3 - Rax_ay_bx_by*mu;

dRb_x_ay_ax = sigma_by*Iby_x_ay_ax*P1 - Rb_x_ay_ax*mu;
dRy_x_ay_ax = sigma_by*Iby_x_ay_ax*P2 - Ry_x_ay_ax*mu;
dRby_x_ay_ax = sigma_by*Iby_x_ay_ax*P3 - Rby_x_ay_ax*mu;

dRb_bx_ay_ax = sigma_by*Iby_bx_ay_ax*P1 - Rb_bx_ay_ax*mu;
dRy_bx_ay_ax = sigma_by*Iby_bx_ay_ax*P2 - Ry_bx_ay_ax*mu;
dRby_bx_ay_ax = sigma_by*Iby_bx_ay_ax*P3 - Rby_bx_ay_ax*mu;

dRa_b_by_ay = sigma_ax*Iax_b_by_ay*P1 - Ra_b_by_ay*mu;
dRx_b_by_ay = sigma_ax*Iax_b_by_ay*P2 - Rx_b_by_ay*mu;
dRax_b_by_ay = sigma_ax*Iax_b_by_ay*P3 - Rax_b_by_ay*mu;

dRa_bx_by_ay = sigma_ax*Iax_bx_by_ay*P1 - Ra_bx_by_ay*mu;
dRx_bx_by_ay = sigma_ax*Iax_bx_by_ay*P2 - Rx_bx_by_ay*mu;
dRax_bx_by_ay = sigma_ax*Iax_bx_by_ay*P3 - Rax_bx_by_ay*mu;

dRa_b_bx_ax = sigma_ay*Iay_b_bx_ax*P1 - Ra_b_bx_ax*mu;
dRy_b_bx_ax = sigma_ay*Iay_b_bx_ax*P2 - Ry_b_bx_ax*mu;
dRay_b_bx_ax = sigma_ay*Iay_b_bx_ax*P3 - Ray_b_bx_ax*mu;

dRa_by_bx_ax = sigma_ay*Iay_by_bx_ax*P1 - Ra_by_bx_ax*mu;
dRy_by_bx_ax = sigma_ay*Iay_by_bx_ax*P2 - Ry_by_bx_ax*mu;
dRay_by_bx_ax = sigma_ay*Iay_by_bx_ax*P3 - Ray_by_bx_ax*mu;

dRb_y_ay_ax = sigma_bx*Ibx_y_ay_ax*P1 - Rb_y_ay_ax*mu;
dRx_y_ay_ax = sigma_bx*Ibx_y_ay_ax*P2 - Rx_y_ay_ax*mu;
dRbx_y_ay_ax = sigma_bx*Ibx_y_ay_ax*P3 - Rbx_y_ay_ax*mu;

dRb_by_ay_ax = sigma_bx*Ibx_by_ay_ax*P1 - Rb_by_ay_ax*mu;
dRx_by_ay_ax = sigma_bx*Ibx_by_ay_ax*P2 - Rx_by_ay_ax*mu;
dRbx_by_ay_ax = sigma_bx*Ibx_by_ay_ax*P3 - Rbx_by_ay_ax*mu;

dRb_by_a_ax_bx  = sigma_by*Iby_a_ax_bx*P1 - Rb_by_a_ax_bx*mu;
dRy_by_a_ax_bx  = sigma_by*Iby_a_ax_bx*P2 - Ry_by_a_ax_bx*mu;
dRby_by_a_ax_bx = sigma_by*Iby_a_ax_bx*P3 - Rby_by_a_ax_bx*mu;

dRb_by_ay_ax_bx  = sigma_by*Iby_ay_ax_bx*P1 - Rb_by_ay_ax_bx*mu;
dRy_by_ay_ax_bx  = sigma_by*Iby_ay_ax_bx*P2 - Ry_by_ay_ax_bx*mu;
dRby_by_ay_ax_bx = sigma_by*Iby_ay_ax_bx*P3 - Rby_by_ay_ax_bx*mu;

dRa_ay_b_ax_bx  = sigma_ay*Iay_b_ax_bx*P1 - Ra_ay_b_ax_bx*mu;
dRy_ay_b_ax_bx  = sigma_ay*Iay_b_ax_bx*P2 - Ry_ay_b_ax_bx*mu;
dRay_ay_b_ax_bx = sigma_ay*Iay_b_ax_bx*P3 - Ray_ay_b_ax_bx*mu;

dRa_ay_by_ax_bx  = sigma_ay*Iay_by_ax_bx*P1 - Ra_ay_by_ax_bx*mu;
dRy_ay_by_ax_bx  = sigma_ay*Iay_by_ax_bx*P2 - Ry_ay_by_ax_bx*mu;
dRay_ay_by_ax_bx = sigma_ay*Iay_by_ax_bx*P3 - Ray_ay_by_ax_bx*mu;

dRb_by_x_ax_ay  = sigma_by*Iby_x_ax_ay*P1 - Rb_by_x_ax_ay*mu;
dRy_by_x_ax_ay  = sigma_by*Iby_x_ax_ay*P2 - Ry_by_x_ax_ay*mu;
dRby_by_x_ax_ay = sigma_by*Iby_x_ax_ay*P3 - Rby_by_x_ax_ay*mu;

dRb_by_bx_ax_ay  = sigma_by*Iby_bx_ax_ay*P1 - Rb_by_bx_ax_ay*mu;
dRy_by_bx_ax_ay  = sigma_by*Iby_bx_ax_ay*P2 - Ry_by_bx_ax_ay*mu;
dRby_by_bx_ax_ay = sigma_by*Iby_bx_ax_ay*P3 - Rby_by_bx_ax_ay*mu;

dRb_bx_y_ax_ay  = sigma_bx*Ibx_y_ax_ay*P1 - Rb_bx_y_ax_ay*mu;
dRx_bx_y_ax_ay  = sigma_bx*Ibx_y_ax_ay*P2 - Rx_bx_y_ax_ay*mu;
dRbx_bx_y_ax_ay = sigma_bx*Ibx_y_ax_ay*P3 - Rbx_bx_y_ax_ay*mu;

dRb_bx_by_ax_ay  = sigma_bx*Ibx_by_ax_ay*P1 - Rb_bx_by_ax_ay*mu;
dRx_bx_by_ax_ay  = sigma_bx*Ibx_by_ax_ay*P2 - Rx_bx_by_ax_ay*mu;
dRbx_bx_by_ax_ay = sigma_bx*Ibx_by_ax_ay*P3 - Rbx_bx_by_ax_ay*mu;

dRa_ay_x_by_bx  = sigma_ay*Iay_x_by_bx*P1 - Ra_ay_x_by_bx*mu;
dRy_ay_x_by_bx  = sigma_ay*Iay_x_by_bx*P2 - Ry_ay_x_by_bx*mu;
dRay_ay_x_by_bx = sigma_ay*Iay_x_by_bx*P3 - Ray_ay_x_by_bx*mu;

dRa_ay_ax_by_bx  = sigma_ay*Iay_ax_by_bx*P1 - Ra_ay_ax_by_bx*mu;
dRy_ay_ax_by_bx  = sigma_ay*Iay_ax_by_bx*P2 - Ry_ay_ax_by_bx*mu;
dRay_ay_ax_by_bx = sigma_ay*Iay_ax_by_bx*P3 - Ray_ay_ax_by_bx*mu;

dRa_ax_y_by_bx  = sigma_ax*Iax_y_by_bx*P1 - Ra_ax_y_by_bx*mu;
dRx_ax_y_by_bx  = sigma_ax*Iax_y_by_bx*P2 - Rx_ax_y_by_bx*mu;
dRax_ax_y_by_bx = sigma_ax*Iax_y_by_bx*P3 - Rax_ax_y_by_bx*mu;

dRa_ax_ay_by_bx  = sigma_ax*Iax_ay_by_bx*P1 - Ra_ax_ay_by_bx*mu;
dRx_ax_ay_by_bx  = sigma_ax*Iax_ay_by_bx*P2 - Rx_ax_ay_by_bx*mu;
dRax_ax_ay_by_bx = sigma_ax*Iax_ay_by_bx*P3 - Rax_ax_ay_by_bx*mu;

dRb_bx_a_ay_by  = sigma_bx*Ibx_a_ay_by*P1 - Rb_bx_a_ay_by*mu;
dRx_bx_a_ay_by  = sigma_bx*Ibx_a_ay_by*P2 - Rx_bx_a_ay_by*mu;
dRbx_bx_a_ay_by = sigma_bx*Ibx_a_ay_by*P3 - Rbx_bx_a_ay_by*mu;

dRb_bx_ax_ay_by  = sigma_bx*Ibx_ax_ay_by*P1 - Rb_bx_ax_ay_by*mu;
dRx_bx_ax_ay_by  = sigma_bx*Ibx_ax_ay_by*P2 - Rx_bx_ax_ay_by*mu;
dRbx_bx_ax_ay_by = sigma_bx*Ibx_ax_ay_by*P3 - Rbx_bx_ax_ay_by*mu;

dRa_ax_b_ay_by  = sigma_ax*Iax_b_ay_by*P1 - Ra_ax_b_ay_by*mu;
dRx_ax_b_ay_by  = sigma_ax*Iax_b_ay_by*P2 - Rx_ax_b_ay_by*mu;
dRax_ax_b_ay_by = sigma_ax*Iax_b_ay_by*P3 - Rax_ax_b_ay_by*mu;

dRa_ax_bx_ay_by  = sigma_ax*Iax_bx_ay_by*P1 - Ra_ax_bx_ay_by*mu;
dRx_ax_bx_ay_by  = sigma_ax*Iax_bx_ay_by*P2 - Rx_ax_bx_ay_by*mu;
dRax_ax_bx_ay_by = sigma_ax*Iax_bx_ay_by*P3 - Rax_ax_bx_ay_by*mu;

dRx_b_by = sigma_ax*Iax_b_by*P2 - Rx_b_by*(lambda_ay*(1 - gamma) + mu);
dRax_b_by = sigma_ax*Iax_b_by*P3 - Rax_b_by*(lambda_ay*(1 - gamma) + mu);

dRa_y_by = sigma_ax*Iax_y_by*P1 - Ra_y_by*(lambda_bx*(1 - gamma) + mu);

dRax_y_by = sigma_ax*Iax_y_by*P3 - Rax_y_by*(lambda_bx*(1 - gamma) + mu);

dRy_b_bx = sigma_ay*Iay_b_bx*P2 - Ry_b_bx*(lambda_ax*(1 - gamma) + mu);
dRay_b_bx = sigma_ay*Iay_b_bx*P3 - Ray_b_bx*(lambda_ax*(1 - gamma) + mu);

dRa_x_bx = sigma_ay*Iay_x_bx*P1 - Ra_x_bx*(lambda_by*(1 - gamma) + mu);

dRay_x_bx = sigma_ay*Iay_x_bx*P3 - Ray_x_bx*(lambda_by*(1 - gamma) + mu);

dRx_a_ay = sigma_bx*Ibx_a_ay*P2 - Rx_a_ay*(lambda_by*(1 - gamma) + mu);
dRbx_a_ay = sigma_bx*Ibx_a_ay*P3 - Rbx_a_ay*(lambda_by*(1 - gamma) + mu);

dRb_y_ay = sigma_bx*Ibx_y_ay*P1 - Rb_y_ay*(lambda_ax*(1 - gamma) + mu); 

dRbx_y_ay = sigma_bx*Ibx_y_ay*P3 - Rbx_y_ay*(lambda_ax*(1 - gamma) + mu);

dRy_a_ax = sigma_by*Iby_a_ax*P2 - Ry_a_ax*(lambda_bx*(1 - gamma) + mu);
dRby_a_ax = sigma_by*Iby_a_ax*P3 - Rby_a_ax*(lambda_bx*(1 - gamma) + mu);

dRb_x_ax = sigma_by*Iby_x_ax*P1 - Rb_x_ax*(lambda_ay*(1 - gamma) + mu);

dRby_x_ax = sigma_by*Iby_x_ax*P3 - Rby_x_ax*(lambda_ay*(1 - gamma) + mu);

dIay_x_b_by = lambda_ay*(1 - gamma)*Rx_b_by - Iay_x_b_by*(sigma_ay + mu);
dIay_ax_b_by = lambda_ay*(1 - gamma)*Rax_b_by - Iay_ax_b_by*(sigma_ay + mu);

dIbx_a_y_by = lambda_bx*(1 - gamma)*Ra_y_by - Ibx_a_y_by*(sigma_bx + mu);
dIbx_ax_y_by = lambda_bx*(1 - gamma)*Rax_y_by - Ibx_ax_y_by*(sigma_bx + mu);

dIax_y_b_bx = lambda_ax*(1 - gamma)*Ry_b_bx - Iax_y_b_bx*(sigma_ax + mu);
dIax_ay_b_bx = lambda_ax*(1 - gamma)*Ray_b_bx - Iax_ay_b_bx*(sigma_ax + mu);

dIby_a_x_bx = lambda_by*(1 - gamma)*Ra_x_bx - Iby_a_x_bx*(sigma_by + mu);
dIby_ay_x_bx = lambda_by*(1 - gamma)*Ray_x_bx - Iby_ay_x_bx*(sigma_by + mu);

dIby_x_a_ay = lambda_by*(1 - gamma)*Rx_a_ay - Iby_x_a_ay*(sigma_by + mu);
dIby_bx_a_ay = lambda_by*(1 - gamma)*Rbx_a_ay - Iby_bx_a_ay*(sigma_by + mu);

dIax_b_y_ay = lambda_ax*(1 - gamma)*Rb_y_ay - Iax_b_y_ay*(sigma_ax + mu);
dIax_bx_y_ay = lambda_ax*(1 - gamma)*Rbx_y_ay - Iax_bx_y_ay*(sigma_ax + mu);

dIbx_y_a_ax = lambda_bx*(1 - gamma)*Ry_a_ax - Ibx_y_a_ax*(sigma_bx + mu);
dIbx_by_a_ax = lambda_bx*(1 - gamma)*Rby_a_ax - Ibx_by_a_ax*(sigma_bx + mu);

dIay_b_x_ax = lambda_ay*(1 - gamma)*Rb_x_ax - Iay_b_x_ax*(sigma_ay + mu);
dIay_by_x_ax = lambda_ay*(1 - gamma)*Rby_x_ax - Iay_by_x_ax*(sigma_ay + mu);

dRa_x_b_by = sigma_ay*Iay_x_b_by*P1 - Ra_x_b_by*mu;
dRy_x_b_by = sigma_ay*Iay_x_b_by*P2 - Ry_x_b_by*mu;
dRay_x_b_by = sigma_ay*Iay_x_b_by*P3 - Ray_x_b_by*mu;

dRa_ax_b_by = sigma_ay*Iay_ax_b_by*P1 - Ra_ax_b_by*mu;
dRy_ax_b_by = sigma_ay*Iay_ax_b_by*P2 - Ry_ax_b_by*mu;
dRay_ax_b_by = sigma_ay*Iay_ax_b_by*P3 - Ray_ax_b_by*mu;

dRb_a_y_by = sigma_bx*Ibx_a_y_by*P1 - Rb_a_y_by*mu;
dRx_a_y_by = sigma_bx*Ibx_a_y_by*P2 - Rx_a_y_by*mu;
dRbx_a_y_by = sigma_bx*Ibx_a_y_by*P3 - Rbx_a_y_by*mu;

dRb_ax_y_by = sigma_bx*Ibx_ax_y_by*P1 - Rb_ax_y_by*mu;
dRx_ax_y_by = sigma_bx*Ibx_ax_y_by*P2 - Rx_ax_y_by*mu;
dRbx_ax_y_by = sigma_bx*Ibx_ax_y_by*P3 - Rbx_ax_y_by*mu;

dRa_y_b_bx = sigma_ax*Iax_y_b_bx*P1 - Ra_y_b_bx*mu;
dRx_y_b_bx = sigma_ax*Iax_y_b_bx*P2 - Rx_y_b_bx*mu;
dRax_y_b_bx = sigma_ax*Iax_y_b_bx*P3 - Rax_y_b_bx*mu;

dRa_ay_b_bx = sigma_ax*Iax_ay_b_bx*P1 - Ra_ay_b_bx*mu;
dRx_ay_b_bx = sigma_ax*Iax_ay_b_bx*P2 - Rx_ay_b_bx*mu;
dRax_ay_b_bx = sigma_ax*Iax_ay_b_bx*P3 - Rax_ay_b_bx*mu;

dRb_a_x_bx = sigma_by*Iby_a_x_bx*P1 - Rb_a_x_bx*mu;
dRy_a_x_bx = sigma_by*Iby_a_x_bx*P2 - Ry_a_x_bx*mu;
dRby_a_x_bx = sigma_by*Iby_a_x_bx*P3 - Rby_a_x_bx*mu;

dRb_ay_x_bx = sigma_by*Iby_ay_x_bx*P1 - Rb_ay_x_bx*mu;
dRy_ay_x_bx = sigma_by*Iby_ay_x_bx*P2 - Ry_ay_x_bx*mu;
dRby_ay_x_bx = sigma_by*Iby_ay_x_bx*P3 - Rby_ay_x_bx*mu;

dRb_x_a_ay = sigma_by*Iby_x_a_ay*P1 - Rb_x_a_ay*mu;
dRy_x_a_ay = sigma_by*Iby_x_a_ay*P2 - Ry_x_a_ay*mu;
dRby_x_a_ay = sigma_by*Iby_x_a_ay*P3 - Rby_x_a_ay*mu;

dRb_bx_a_ay = sigma_by*Iby_bx_a_ay*P1 - Rb_bx_a_ay*mu;
dRy_bx_a_ay = sigma_by*Iby_bx_a_ay*P2 - Ry_bx_a_ay*mu;
dRby_bx_a_ay = sigma_by*Iby_bx_a_ay*P3 - Rby_bx_a_ay*mu;

dRa_b_y_ay = sigma_ax*Iax_b_y_ay*P1 - Ra_b_y_ay*mu;
dRx_b_y_ay = sigma_ax*Iax_b_y_ay*P2 - Rx_b_y_ay*mu;
dRax_b_y_ay = sigma_ax*Iax_b_y_ay*P3 - Rax_b_y_ay*mu;

dRa_bx_y_ay = sigma_ax*Iax_bx_y_ay*P1 - Ra_bx_y_ay*mu;
dRx_bx_y_ay = sigma_ax*Iax_bx_y_ay*P2 - Rx_bx_y_ay*mu;
dRax_bx_y_ay = sigma_ax*Iax_bx_y_ay*P3 - Rax_bx_y_ay*mu;

dRb_y_a_ax = sigma_bx*Ibx_y_a_ax*P1 - Rb_y_a_ax*mu;
dRx_y_a_ax = sigma_bx*Ibx_y_a_ax*P2 - Rx_y_a_ax*mu;
dRbx_y_a_ax = sigma_bx*Ibx_y_a_ax*P3 - Rbx_y_a_ax*mu;

dRb_by_a_ax = sigma_bx*Ibx_by_a_ax*P1 - Rb_by_a_ax*mu;
dRx_by_a_ax = sigma_bx*Ibx_by_a_ax*P2 - Rx_by_a_ax*mu;
dRbx_by_a_ax = sigma_bx*Ibx_by_a_ax*P3 - Rbx_by_a_ax*mu;

dRa_b_x_ax = sigma_ay*Iay_b_x_ax*P1 - Ra_b_x_ax*mu;
dRy_b_x_ax = sigma_ay*Iay_b_x_ax*P2 - Ry_b_x_ax*mu;
dRay_b_x_ax = sigma_ay*Iay_b_x_ax*P3 - Ray_b_x_ax*mu;

dRa_by_x_ax = sigma_ay*Iay_by_x_ax*P1 - Ra_by_x_ax*mu;
dRy_by_x_ax = sigma_ay*Iay_by_x_ax*P2 - Ry_by_x_ax*mu;
dRay_by_x_ax = sigma_ay*Iay_by_x_ax*P3 - Ray_by_x_ax*mu;

dIax_y_ay = lambda_ax*Ry_ay*(1 - gamma) - Iax_y_ay*(sigma_ax + mu);
dIax_b_bx = lambda_ax*Rb_bx*(1 - gamma) - Iax_b_bx*(sigma_ax + mu);

dIay_x_ax = lambda_ay*Rx_ax*(1 - gamma) - Iay_x_ax*(sigma_ay + mu);
dIay_b_by = lambda_ay*Rb_by*(1 - gamma) - Iay_b_by*(sigma_ay + mu);

dIbx_a_ax = lambda_bx*Ra_ax*(1 - gamma) - Ibx_a_ax*(sigma_bx + mu);
dIbx_y_by = lambda_bx*Ry_by*(1 - gamma) - Ibx_y_by*(sigma_bx + mu);

dIby_a_ay = lambda_by*Ra_ay*(1 - gamma) - Iby_a_ay*(sigma_by + mu);
dIby_x_bx = lambda_by*Rx_bx*(1 - gamma) - Iby_x_bx*(sigma_by + mu);

dRa_x_ax = sigma_ay*Iay_x_ax*P1 - Ra_x_ax*(lambda_by + mu);
dRa_y_ay = sigma_ax*Iax_y_ay*P1 - Ra_y_ay*(lambda_bx + mu);
dRb_x_bx = sigma_by*Iby_x_bx*P1 - Rb_x_bx*(lambda_ay + mu);
dRb_y_by = sigma_bx*Ibx_y_by*P1 - Rb_y_by*(lambda_ax + mu);

dRx_a_ax = sigma_bx*Ibx_a_ax*P2 - Rx_a_ax*(lambda_by + mu);
dRy_a_ay = sigma_by*Iby_a_ay*P2 - Ry_a_ay*(lambda_bx + mu);
dRx_b_bx = sigma_ax*Iax_b_bx*P2 - Rx_b_bx*(lambda_ay + mu);
dRy_b_by = sigma_ay*Iay_b_by*P2 - Ry_b_by*(lambda_ax + mu);

dRa_b_bx = sigma_ax*Iax_b_bx*P1 + sigma_ay*Iay_b_bx*P1 - Ra_b_bx*mu;
dRa_b_by = sigma_ay*Iay_b_by*P1 + sigma_ax*Iax_b_by*P1 - Ra_b_by*mu; 
dRb_a_ax = sigma_bx*Ibx_a_ax*P1 + sigma_by*Iby_a_ax*P1 - Rb_a_ax*mu;
dRb_a_ay = sigma_by*Iby_a_ay*P1 + sigma_bx*Ibx_a_ay*P1 - Rb_a_ay*mu;
dRx_y_ay = sigma_ax*Iax_y_ay*P2 + sigma_bx*Ibx_y_ay*P2 - Rx_y_ay*mu;
dRy_x_ax = sigma_ay*Iay_x_ax*P2 + sigma_by*Iby_x_ax*P2 - Ry_x_ax*mu;
dRx_y_by = sigma_bx*Ibx_y_by*P2 + sigma_ax*Iax_y_by*P2 - Rx_y_by*mu;
dRy_x_bx = sigma_by*Iby_x_bx*P2 + sigma_ay*Iay_x_bx*P2 - Ry_x_bx*mu;

dRbx_a_ax = sigma_bx*Ibx_a_ax*P3 - Rbx_a_ax*(lambda_by*(1 - gamma) + mu);
dRay_x_ax = sigma_ay*Iay_x_ax*P3 - Ray_x_ax*(lambda_by*(1 - gamma) + mu);
dRby_a_ay = sigma_by*Iby_a_ay*P3 - Rby_a_ay*(lambda_bx*(1 - gamma) + mu);
dRax_y_ay = sigma_ax*Iax_y_ay*P3 - Rax_y_ay*(lambda_bx*(1 - gamma) + mu);
dRby_x_bx = sigma_by*Iby_x_bx*P3 - Rby_x_bx*(lambda_ay*(1 - gamma) + mu);
dRax_b_bx = sigma_ax*Iax_b_bx*P3 - Rax_b_bx*(lambda_ay*(1 - gamma) + mu);
dRbx_y_by = sigma_bx*Ibx_y_by*P3 - Rbx_y_by*(lambda_ax*(1 - gamma) + mu);
dRay_b_by = sigma_ay*Iay_b_by*P3 - Ray_b_by*(lambda_ax*(1 - gamma) + mu);

dIax_b_y_by = lambda_ax*Rb_y_by - Iax_b_y_by*(sigma_ax + mu);
dIay_b_x_bx = lambda_ay*Rb_x_bx - Iay_b_x_bx*(sigma_ay + mu);
dIbx_a_y_ay = lambda_bx*Ra_y_ay - Ibx_a_y_ay*(sigma_bx + mu);
dIby_a_x_ax = lambda_by*Ra_x_ax - Iby_a_x_ax*(sigma_by + mu);

dIby_x_a_ax = lambda_by*Rx_a_ax - (sigma_by + mu)*Iby_x_a_ax;
dIbx_y_a_ay = lambda_bx*Ry_a_ay - (sigma_bx + mu)*Ibx_y_a_ay;
dIay_x_b_bx = lambda_ay*Rx_b_bx - (sigma_ay + mu)*Iay_x_b_bx;
dIax_y_b_by = lambda_ax*Ry_b_by - (sigma_ax + mu)*Iax_y_b_by;

dIax_ay_b_by = lambda_ax*(1 - gamma)*Ray_b_by - Iax_ay_b_by*(sigma_ax + mu);
dIax_bx_y_by = lambda_ax*(1 - gamma)*Rbx_y_by - Iax_bx_y_by*(sigma_ax + mu);

dIay_ax_b_bx = lambda_ay*(1 - gamma)*Rax_b_bx - Iay_ax_b_bx*(sigma_ay + mu);
dIay_by_x_bx = lambda_ay*(1 - gamma)*Rby_x_bx - Iay_by_x_bx*(sigma_ay + mu);

dIbx_ax_y_ay = lambda_bx*(1 - gamma)*Rax_y_ay - Ibx_ax_y_ay*(sigma_bx + mu);
dIbx_by_a_ay = lambda_bx*(1 - gamma)*Rby_a_ay - Ibx_by_a_ay*(sigma_bx + mu);

dIby_ay_x_ax = lambda_by*(1 - gamma)*Ray_x_ax - Iby_ay_x_ax*(sigma_by + mu);
dIby_bx_a_ax = lambda_by*(1 - gamma)*Rbx_a_ax - Iby_bx_a_ax*(sigma_by + mu);

dRa_b_y_by   = sigma_ax*Iax_b_y_by*P1 - Ra_b_y_by*mu;
dRx_b_y_by   = sigma_ax*Iax_b_y_by*P2 - Rx_b_y_by*mu;
dRax_b_y_by  = sigma_ax*Iax_b_y_by*P3 - Rax_b_y_by*mu;

dRa_b_x_bx   = sigma_ay*Iay_b_x_bx*P1 - Ra_b_x_bx*mu;
dRy_b_x_bx   = sigma_ay*Iay_b_x_bx*P2 - Ry_b_x_bx*mu;
dRay_b_x_bx  = sigma_ay*Iay_b_x_bx*P3 - Ray_b_x_bx*mu;

dRb_a_y_ay   = sigma_bx*Ibx_a_y_ay*P1 - Rb_a_y_ay*mu;
dRx_a_y_ay   = sigma_bx*Ibx_a_y_ay*P2 - Rx_a_y_ay*mu;
dRbx_a_y_ay  = sigma_bx*Ibx_a_y_ay*P3 - Rbx_a_y_ay*mu;

dRb_a_x_ax   = sigma_by*Iby_a_x_ax*P1 - Rb_a_x_ax*mu;
dRy_a_x_ax   = sigma_by*Iby_a_x_ax*P2 - Ry_a_x_ax*mu;
dRby_a_x_ax  = sigma_by*Iby_a_x_ax*P3 - Rby_a_x_ax*mu;

dRb_x_a_ax   = sigma_by*Iby_x_a_ax*P1 - Rb_x_a_ax*mu;
dRy_x_a_ax   = sigma_by*Iby_x_a_ax*P2 - Ry_x_a_ax*mu;
dRby_x_a_ax  = sigma_by*Iby_x_a_ax*P3 - Rby_x_a_ax*mu;

dRb_y_a_ay   = sigma_bx*Ibx_y_a_ay*P1 - Rb_y_a_ay*mu;
dRx_y_a_ay   = sigma_bx*Ibx_y_a_ay*P2 - Rx_y_a_ay*mu;
dRbx_y_a_ay  = sigma_bx*Ibx_y_a_ay*P3 - Rbx_y_a_ay*mu;

dRa_x_b_bx   = sigma_ay*Iay_x_b_bx*P1 - Ra_x_b_bx*mu;
dRy_x_b_bx   = sigma_ay*Iay_x_b_bx*P2 - Ry_x_b_bx*mu;
dRay_x_b_bx  = sigma_ay*Iay_x_b_bx*P3 - Ray_x_b_bx*mu;

dRa_y_b_by   = sigma_ax*Iax_y_b_by*P1 - Ra_y_b_by*mu;
dRx_y_b_by   = sigma_ax*Iax_y_b_by*P2 - Rx_y_b_by*mu;
dRax_y_b_by  = sigma_ax*Iax_y_b_by*P3 - Rax_y_b_by*mu;

dRa_ay_b_by = sigma_ax*Iax_ay_b_by*P1 - Ra_ay_b_by*mu;
dRx_ay_b_by = sigma_ax*Iax_ay_b_by*P2 - Rx_ay_b_by*mu;
dRax_ay_b_by = sigma_ax*Iax_ay_b_by*P3 - Rax_ay_b_by*mu;

dRa_bx_y_by = sigma_ax*Iax_bx_y_by*P1 - Ra_bx_y_by*mu;
dRx_bx_y_by = sigma_ax*Iax_bx_y_by*P2 - Rx_bx_y_by*mu;
dRax_bx_y_by = sigma_ax*Iax_bx_y_by*P3 - Rax_bx_y_by*mu;

dRa_ax_b_bx = sigma_ay*Iay_ax_b_bx*P1 - Ra_ax_b_bx*mu;
dRy_ax_b_bx = sigma_ay*Iay_ax_b_bx*P2 - Ry_ax_b_bx*mu;
dRay_ax_b_bx = sigma_ay*Iay_ax_b_bx*P3 - Ray_ax_b_bx*mu;

dRa_by_x_bx = sigma_ay*Iay_by_x_bx*P1 - Ra_by_x_bx*mu;
dRy_by_x_bx = sigma_ay*Iay_by_x_bx*P2 - Ry_by_x_bx*mu;
dRay_by_x_bx = sigma_ay*Iay_by_x_bx*P3 - Ray_by_x_bx*mu;

dRb_ax_y_ay = sigma_bx*Ibx_ax_y_ay*P1 - Rb_ax_y_ay*mu;
dRx_ax_y_ay = sigma_bx*Ibx_ax_y_ay*P2 - Rx_ax_y_ay*mu;
dRbx_ax_y_ay = sigma_bx*Ibx_ax_y_ay*P3 - Rbx_ax_y_ay*mu;

dRb_by_a_ay = sigma_bx*Ibx_by_a_ay*P1 - Rb_by_a_ay*mu;
dRx_by_a_ay = sigma_bx*Ibx_by_a_ay*P2 - Rx_by_a_ay*mu;
dRbx_by_a_ay = sigma_bx*Ibx_by_a_ay*P3 - Rbx_by_a_ay*mu;

dRb_ay_x_ax = sigma_by*Iby_ay_x_ax*P1 - Rb_ay_x_ax*mu;
dRy_ay_x_ax = sigma_by*Iby_ay_x_ax*P2 - Ry_ay_x_ax*mu;
dRby_ay_x_ax = sigma_by*Iby_ay_x_ax*P3 - Rby_ay_x_ax*mu;

dRb_bx_a_ax = sigma_by*Iby_bx_a_ax*P1 - Rb_bx_a_ax*mu;
dRy_bx_a_ax = sigma_by*Iby_bx_a_ax*P2 - Ry_bx_a_ax*mu;
dRby_bx_a_ax = sigma_by*Iby_bx_a_ax*P3 - Rby_bx_a_ax*mu;


dydt = [
    dS;
    dIax;
    dIay;
    dIbx;
    dIby;
    dRa;
    dRb;
    dRx;
    dRy;
    dRax_P3;
    dRay_P3;
    dRbx_P3;
    dRby_P3;
    dIax_y;
    dIax_b;
    dIay_x;
    dIay_b;
    dIbx_y;
    dIbx_a;
    dIby_x;
    dIby_a;
    dRa_x;
    dRa_y;
    dRb_x;
    dRb_y;
    dRx_y;
    dRa_b;
    dRx_a;
    dRy_a;
    dRx_b;
    dRy_b;
    dRax_y;
    dRax_b;
    dRay_x;
    dRay_b;
    dRbx_y;
    dRbx_a;
    dRby_x;
    dRby_a;
    dIby_a_x;
    dIbx_a_y;
    dIay_b_x;
    dIax_b_y;
    dIby_x_a;
    dIbx_y_a;
    dIay_x_b;
    dIax_y_b;
    dIax_ay_b;
    dIax_bx_y;
    dIay_ax_b;
    dIay_by_x;
    dIbx_ax_y;
    dIbx_by_a;
    dIby_ay_x;
    dIby_bx_a;
    dRa_b_y;
    dRx_b_y;
    dRax_b_y;
    dRa_b_x;
    dRy_b_x;
    dRay_b_x;
    dRb_a_y;
    dRx_a_y;
    dRbx_a_y;
    dRb_a_x;
    dRy_a_x;
    dRby_a_x;
    dRa_y_b;
    dRx_y_b;
    dRax_y_b;
    dRa_x_b;
    dRy_x_b;
    dRay_x_b;
    dRb_y_a;
    dRx_y_a;
    dRbx_y_a;
    dRb_x_a;
    dRy_x_a;
    dRby_x_a;
    dRa_ay_b;
    dRx_ay_b;
    dRax_ay_b;
    dRa_bx_y;
    dRx_bx_y;
    dRax_bx_y;
    dRa_ax_b;
    dRy_ax_b;
    dRay_ax_b;
    dRa_by_x;
    dRy_by_x;
    dRay_by_x;
    dRb_ax_y;
    dRx_ax_y;
    dRbx_ax_y;
    dRb_by_a;
    dRx_by_a;
    dRbx_by_a;
    dRb_ay_x;
    dRy_ay_x;
    dRby_ay_x;
    dRb_bx_a;
    dRy_bx_a;
    dRby_bx_a;
    dIay_ax;
    dIbx_ax;
    dIax_ay;
    dIby_ay;
    dIax_bx;
    dIby_bx;
    dIay_by;
    dIbx_by;
    dIax_by;
    dIay_bx;
    dIbx_ay;
    dIby_ax;
    dRa_ax;
    dRy_ax;
    dRb_ax;
    dRx_ax;
    dRa_ay;
    dRx_ay;
    dRb_ay;
    dRy_ay;
    dRa_bx;
    dRx_bx;
    dRb_bx;
    dRy_bx;
    dRa_by;
    dRy_by;
    dRb_by;
    dRx_by;
    dRbx_ax;
    dRay_ax;
    dRbx_by;
    dRby_ay;
    dRax_bx;
    dRax_ay;
    dRby_bx;
    dRay_by;
    dRax_by;
    dRay_bx;
    dRby_ax;
    dRbx_ay;
    dIax_ay_bx;
    dIay_ax_by;
    dIbx_ax_by;
    dIby_ay_bx;
    dIay_by_ax;
    dIbx_by_ax;
    dIax_bx_ay;
    dIby_bx_ay;
    dRa_ay_bx;
    dRx_ay_bx;
    dRax_ay_bx;
    dRa_ax_by;
    dRy_ax_by;
    dRay_ax_by;
    dRb_ax_by;
    dRx_ax_by;
    dRbx_ax_by;
    dRb_ay_bx;
    dRy_ay_bx;
    dRby_ay_bx;
    dRa_by_ax;
    dRy_by_ax;
    dRay_by_ax;
    dRb_by_ax;
    dRx_by_ax;
    dRbx_by_ax;
    dRa_bx_ay;
    dRx_bx_ay;
    dRax_bx_ay;
    dRb_bx_ay;
    dRy_bx_ay;
    dRby_bx_ay;
    dIax_b_ay_bx;
    dIax_y_ay_bx;
    dIax_by_ay_bx;
    dIay_b_ax_by;
    dIay_x_ax_by;
    dIay_bx_ax_by;
    dIbx_a_ax_by;
    dIbx_y_ax_by;
    dIbx_ay_ax_by;
    dIby_a_ay_bx;
    dIby_x_ay_bx;
    dIby_ax_ay_bx;
    dIay_b_by_ax;
    dIay_x_by_ax;
    dIay_bx_by_ax;
    dIbx_a_by_ax;
    dIbx_y_by_ax;
    dIbx_ay_by_ax;
    dIax_b_bx_ay;
    dIax_y_bx_ay;
    dIax_by_bx_ay;
    dIby_a_bx_ay;
    dIby_x_bx_ay;
    dIby_ax_bx_ay;
    dRa_b_ay_bx;
    dRx_b_ay_bx;
    dRax_b_ay_bx;
    dRa_y_ay_bx;
    dRx_y_ay_bx;
    dRax_y_ay_bx;
    dRa_by_ay_bx;
    dRx_by_ay_bx;
    dRax_by_ay_bx;
    dRa_b_ax_by;
    dRy_b_ax_by;
    dRay_b_ax_by;
    dRa_x_ax_by;
    dRy_x_ax_by;
    dRay_x_ax_by;
    dRa_bx_ax_by;
    dRy_bx_ax_by;
    dRay_bx_ax_by;
    dRb_a_ax_by;
    dRx_a_ax_by;
    dRbx_a_ax_by;
    dRb_y_ax_by;
    dRx_y_ax_by;
    dRbx_y_ax_by;
    dRb_ay_ax_by;
    dRx_ay_ax_by;
    dRbx_ay_ax_by;
    dRb_a_ay_bx;
    dRy_a_ay_bx;
    dRby_a_ay_bx;
    dRb_x_ay_bx;
    dRy_x_ay_bx;
    dRby_x_ay_bx;
    dRb_ax_ay_bx;
    dRy_ax_ay_bx;
    dRby_ax_ay_bx;
    dRa_b_bx_ay;
    dRx_b_bx_ay;
    dRax_b_bx_ay;
    dRa_y_bx_ay;
    dRx_y_bx_ay;
    dRax_y_bx_ay;
    dRa_by_bx_ay;
    dRx_by_bx_ay;
    dRax_by_bx_ay;
    dRa_b_by_ax;
    dRy_b_by_ax;
    dRay_b_by_ax;
    dRa_x_by_ax;
    dRy_x_by_ax;
    dRay_x_by_ax;
    dRa_bx_by_ax;
    dRy_bx_by_ax;
    dRay_bx_by_ax;
    dRb_a_by_ax;
    dRx_a_by_ax;
    dRbx_a_by_ax;
    dRb_y_by_ax;
    dRx_y_by_ax;
    dRbx_y_by_ax;
    dRb_ay_by_ax;
    dRx_ay_by_ax;
    dRbx_ay_by_ax;
    dRb_a_bx_ay;
    dRy_a_bx_ay;
    dRby_a_bx_ay;
    dRb_x_bx_ay;
    dRy_x_bx_ay;
    dRby_x_bx_ay;
    dRb_ax_bx_ay;
    dRy_ax_bx_ay;
    dRby_ax_bx_ay;
    dIax_b_ay;
    dIax_by_ay;
    dIax_y_bx;
    dIax_bx_by;
    dIay_b_ax;
    dIay_bx_ax;
    dIay_x_by;
    dIay_bx_by;
    dIbx_y_ax;
    dIbx_ay_ax;
    dIbx_a_by;
    dIbx_by_ay;
    dIay_ax_bx;
    dIby_ax_bx;
    dIbx_ax_ay;
    dIby_ax_ay;
    dIax_by_bx;
    dIay_by_bx;
    dIax_ay_by;
    dIbx_ay_by;
    dIby_a_bx;
    dIby_bx_ax;
    dIby_x_ay;
    dIby_ay_ax;
    dIax_b_by;
    dIax_y_by;
    dIay_b_bx;
    dIay_x_bx;
    dIbx_a_ay;
    dIbx_y_ay;
    dIby_a_ax;
    dIby_x_ax;
    dRa_b_ay;
    dRx_b_ay;
    dRax_b_ay;
    dRa_by_ay;
    dRx_by_ay;
    dRax_by_ay;
    dRa_y_bx;
    dRx_y_bx;
    dRax_y_bx;
    dRa_bx_by;
    dRx_bx_by;
    dRax_bx_by;
    dRa_b_ax;
    dRy_b_ax;
    dRay_b_ax;
    dRa_bx_ax;
    dRy_bx_ax;
    dRay_bx_ax;
    dRa_x_by;
    dRy_x_by;
    dRay_x_by;
    dRy_bx_by;
    dRay_bx_by;
    dRb_y_ax;
    dRx_y_ax;
    dRbx_y_ax;
    dRb_ay_ax;
    dRx_ay_ax;
    dRbx_ay_ax;
    dRb_a_by;
    dRx_a_by;
    dRbx_a_by;
    dRb_by_ay;
    dRbx_by_ay;
    dRb_a_bx;
    dRy_a_bx;
    dRby_a_bx;
    dRb_bx_ax;
    dRby_bx_ax;
    dRb_x_ay;
    dRy_x_ay;
    dRby_x_ay;
    dRy_ay_ax;
    dRby_ay_ax;
    dRa_ax_bx;
    dRy_ax_bx;
    dRay_ax_bx;
    dRb_ax_bx;
    dRby_ax_bx;
    dRb_ax_ay;
    dRx_ax_ay;
    dRbx_ax_ay;
    dRy_ax_ay;
    dRby_ax_ay;
    dRa_by_bx;
    dRx_by_bx;
    dRax_by_bx;
    dRy_by_bx;
    dRay_by_bx;
    dRa_ay_by;
    dRx_ay_by;
    dRax_ay_by;
    dRb_ay_by;
    dRbx_ay_by;
    dIbx_a_by_ay;
    dIbx_ax_by_ay;
    dIay_x_bx_by;
    dIay_ax_bx_by;
    dIby_a_bx_ax;
    dIby_ay_bx_ax;
    dIax_y_bx_by;
    dIax_ay_bx_by;
    dIby_x_ay_ax;
    dIby_bx_ay_ax;
    dIax_b_by_ay;
    dIax_bx_by_ay;
    dIay_b_bx_ax;
    dIay_by_bx_ax;
    dIbx_y_ay_ax;
    dIbx_by_ay_ax;
    dIby_a_ax_bx;
    dIby_ay_ax_bx;
    dIay_b_ax_bx;
    dIay_by_ax_bx;
    dIby_x_ax_ay;
    dIby_bx_ax_ay;
    dIbx_y_ax_ay;
    dIbx_by_ax_ay;
    dIay_x_by_bx;
    dIay_ax_by_bx;
    dIax_y_by_bx;
    dIax_ay_by_bx;
    dIbx_a_ay_by;
    dIbx_ax_ay_by;
    dIax_b_ay_by;
    dIax_bx_ay_by;
    dRb_a_by_ay;
    dRx_a_by_ay;
    dRbx_a_by_ay;
    dRb_ax_by_ay;
    dRx_ax_by_ay;
    dRbx_ax_by_ay;
    dRa_x_bx_by;
    dRy_x_bx_by;
    dRay_x_bx_by;
    dRa_ax_bx_by;
    dRy_ax_bx_by;
    dRay_ax_bx_by;
    dRb_a_bx_ax;
    dRy_a_bx_ax;
    dRby_a_bx_ax;
    dRb_ay_bx_ax;
    dRy_ay_bx_ax;
    dRby_ay_bx_ax;
    dRa_y_bx_by;
    dRx_y_bx_by;
    dRax_y_bx_by;
    dRa_ay_bx_by;
    dRx_ay_bx_by;
    dRax_ay_bx_by;
    dRb_x_ay_ax;
    dRy_x_ay_ax;
    dRby_x_ay_ax;
    dRb_bx_ay_ax;
    dRy_bx_ay_ax;
    dRby_bx_ay_ax;
    dRa_b_by_ay;
    dRx_b_by_ay;
    dRax_b_by_ay;
    dRa_bx_by_ay;
    dRx_bx_by_ay;
    dRax_bx_by_ay;
    dRa_b_bx_ax;
    dRy_b_bx_ax;
    dRay_b_bx_ax;
    dRa_by_bx_ax;
    dRy_by_bx_ax;
    dRay_by_bx_ax;
    dRb_y_ay_ax;
    dRx_y_ay_ax;
    dRbx_y_ay_ax;
    dRb_by_ay_ax;
    dRx_by_ay_ax;
    dRbx_by_ay_ax;
    dRb_by_a_ax_bx;
    dRy_by_a_ax_bx;
    dRby_by_a_ax_bx;
    dRb_by_ay_ax_bx;
    dRy_by_ay_ax_bx;
    dRby_by_ay_ax_bx;
    dRa_ay_b_ax_bx;
    dRy_ay_b_ax_bx;
    dRay_ay_b_ax_bx;
    dRa_ay_by_ax_bx;
    dRy_ay_by_ax_bx;
    dRay_ay_by_ax_bx;
    dRb_by_x_ax_ay;
    dRy_by_x_ax_ay;
    dRby_by_x_ax_ay;
    dRb_by_bx_ax_ay;
    dRy_by_bx_ax_ay;
    dRby_by_bx_ax_ay;
    dRb_bx_y_ax_ay;
    dRx_bx_y_ax_ay;
    dRbx_bx_y_ax_ay;
    dRb_bx_by_ax_ay;
    dRx_bx_by_ax_ay;
    dRbx_bx_by_ax_ay;
    dRa_ay_x_by_bx;
    dRy_ay_x_by_bx;
    dRay_ay_x_by_bx;
    dRa_ay_ax_by_bx;
    dRy_ay_ax_by_bx;
    dRay_ay_ax_by_bx;
    dRa_ax_y_by_bx;
    dRx_ax_y_by_bx;
    dRax_ax_y_by_bx;
    dRa_ax_ay_by_bx;
    dRx_ax_ay_by_bx;
    dRax_ax_ay_by_bx;
    dRb_bx_a_ay_by;
    dRx_bx_a_ay_by;
    dRbx_bx_a_ay_by;
    dRb_bx_ax_ay_by;
    dRx_bx_ax_ay_by;
    dRbx_bx_ax_ay_by;
    dRa_ax_b_ay_by;
    dRx_ax_b_ay_by;
    dRax_ax_b_ay_by;
    dRa_ax_bx_ay_by;
    dRx_ax_bx_ay_by;
    dRax_ax_bx_ay_by;
    dRa_b_by;
    dRx_b_by;
    dRax_b_by;
    dRa_y_by;
    dRx_y_by;
    dRax_y_by;
    dRa_b_bx;
    dRy_b_bx;
    dRay_b_bx;
    dRa_x_bx;
    dRy_x_bx;
    dRay_x_bx;
    dRb_a_ay;
    dRx_a_ay;
    dRbx_a_ay;
    dRb_y_ay;
    dRx_y_ay;
    dRbx_y_ay;
    dRb_a_ax;
    dRy_a_ax;
    dRby_a_ax;
    dRb_x_ax;
    dRy_x_ax;
    dRby_x_ax;
    dIay_x_b_by;
    dIay_ax_b_by;
    dIbx_a_y_by;
    dIbx_ax_y_by;
    dIax_y_b_bx;
    dIax_ay_b_bx;
    dIby_a_x_bx;
    dIby_ay_x_bx;
    dIby_x_a_ay;
    dIby_bx_a_ay;
    dIax_b_y_ay;
    dIax_bx_y_ay;
    dIbx_y_a_ax;
    dIbx_by_a_ax;
    dIay_b_x_ax;
    dIay_by_x_ax;
    dRa_x_b_by;
    dRy_x_b_by;
    dRay_x_b_by;
    dRa_ax_b_by;
    dRy_ax_b_by;
    dRay_ax_b_by;
    dRb_a_y_by;
    dRx_a_y_by;
    dRbx_a_y_by;
    dRb_ax_y_by;
    dRx_ax_y_by;
    dRbx_ax_y_by;
    dRa_y_b_bx;
    dRx_y_b_bx;
    dRax_y_b_bx;
    dRa_ay_b_bx;
    dRx_ay_b_bx;
    dRax_ay_b_bx;
    dRb_a_x_bx;
    dRy_a_x_bx;
    dRby_a_x_bx;
    dRb_ay_x_bx;
    dRy_ay_x_bx;
    dRby_ay_x_bx;
    dRb_x_a_ay;
    dRy_x_a_ay;
    dRby_x_a_ay;
    dRb_bx_a_ay;
    dRy_bx_a_ay;
    dRby_bx_a_ay;
    dRa_b_y_ay;
    dRx_b_y_ay;
    dRax_b_y_ay;
    dRa_bx_y_ay;
    dRx_bx_y_ay;
    dRax_bx_y_ay;
    dRb_y_a_ax;
    dRx_y_a_ax;
    dRbx_y_a_ax;
    dRb_by_a_ax;
    dRx_by_a_ax;
    dRbx_by_a_ax;
    dRa_b_x_ax;
    dRy_b_x_ax;
    dRay_b_x_ax;
    dRa_by_x_ax;
    dRy_by_x_ax;
    dRay_by_x_ax;
    dIax_y_ay;
    dIax_b_bx;
    dIay_x_ax;
    dIay_b_by;
    dIbx_a_ax;
    dIbx_y_by;
    dIby_a_ay;
    dIby_x_bx;
    dRa_x_ax;
    dRa_y_ay;
    dRb_x_bx;
    dRb_y_by;
    dRx_a_ax;
    dRy_a_ay;
    dRx_b_bx;
    dRy_b_by;
    dRbx_a_ax;
    dRay_x_ax;
    dRby_a_ay;
    dRax_y_ay;
    dRby_x_bx;
    dRax_b_bx;
    dRbx_y_by;
    dRay_b_by;
    dIax_b_y_by;
    dIay_b_x_bx;
    dIbx_a_y_ay;
    dIby_a_x_ax;
    dIby_x_a_ax;
    dIbx_y_a_ay;
    dIay_x_b_bx;
    dIax_y_b_by;
    dIax_ay_b_by;
    dIax_bx_y_by;
    dIay_ax_b_bx;
    dIay_by_x_bx;
    dIbx_ax_y_ay;
    dIbx_by_a_ay;
    dIby_ay_x_ax;
    dIby_bx_a_ax;
    dRa_b_y_by;
    dRx_b_y_by;
    dRax_b_y_by;
    dRa_b_x_bx;
    dRy_b_x_bx;
    dRay_b_x_bx;
    dRb_a_y_ay;
    dRx_a_y_ay;
    dRbx_a_y_ay;
    dRb_a_x_ax;
    dRy_a_x_ax;
    dRby_a_x_ax;
    dRb_x_a_ax;
    dRy_x_a_ax;
    dRby_x_a_ax;
    dRb_y_a_ay;
    dRx_y_a_ay;
    dRbx_y_a_ay;
    dRa_x_b_bx;
    dRy_x_b_bx;
    dRay_x_b_bx;
    dRa_y_b_by;
    dRx_y_b_by;
    dRax_y_b_by;
    dRa_ay_b_by;
    dRx_ay_b_by;
    dRax_ay_b_by;
    dRa_bx_y_by;
    dRx_bx_y_by;
    dRax_bx_y_by;
    dRa_ax_b_bx;
    dRy_ax_b_bx;
    dRay_ax_b_bx;
    dRa_by_x_bx;
    dRy_by_x_bx;
    dRay_by_x_bx;
    dRb_ax_y_ay;
    dRx_ax_y_ay;
    dRbx_ax_y_ay;
    dRb_by_a_ay;
    dRx_by_a_ay;
    dRbx_by_a_ay;
    dRb_ay_x_ax;
    dRy_ay_x_ax;
    dRby_ay_x_ax;
    dRb_bx_a_ax;
    dRy_bx_a_ax;
    dRby_bx_a_ax];
end