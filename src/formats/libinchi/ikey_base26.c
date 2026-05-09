/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*
* The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
*/


 /*
     InChIKey: procedures for base-26 encoding
 */


#ifdef _MSC_VER
#if _MSC_VER > 1000
#pragma warning( disable : 4996 )
#endif
#endif

#include <string.h>
#include <stdio.h>
#include "ikey_base26.h"

#include "mode.h"

#include "bcf_s.h"

 /*
     Triplets

     As the 2^14 (16384) is very close to 26^3 (17576), a triplet of uppercase
     letters A..Z encodes 14 bits with good efficiency.
     For speed, we just tabulate triplets below.

     We should throw away 17576-16384= 1192 triplets.
     These are 676 triplets starting from 'E', the most frequent letter in English
     texts (the other 516 are those started at 'T' , "TAA" to "TTV").
 */

static const char t26[][4] =
{
"AAA","AAB","AAC","AAD","AAE","AAF","AAG","AAH","AAI","AAJ","AAK","AAL","AAM","AAN","AAO","AAP",
"AAQ","AAR","AAS","AAT","AAU","AAV","AAW","AAX","AAY","AAZ","ABA","ABB","ABC","ABD","ABE","ABF",
"ABG","ABH","ABI","ABJ","ABK","ABL","ABM","ABN","ABO","ABP","ABQ","ABR","ABS","ABT","ABU","ABV",
"ABW","ABX","ABY","ABZ","ACA","ACB","ACC","ACD","ACE","ACF","ACG","ACH","ACI","ACJ","ACK","ACL",
"ACM","ACN","ACO","ACP","ACQ","ACR","ACS","ACT","ACU","ACV","ACW","ACX","ACY","ACZ","ADA","ADB",
"ADC","ADD","ADE","ADF","ADG","ADH","ADI","ADJ","ADK","ADL","ADM","ADN","ADO","ADP","ADQ","ADR",
"ADS","ADT","ADU","ADV","ADW","ADX","ADY","ADZ","AEA","AEB","AEC","AED","AEE","AEF","AEG","AEH",
"AEI","AEJ","AEK","AEL","AEM","AEN","AEO","AEP","AEQ","AER","AES","AET","AEU","AEV","AEW","AEX",
"AEY","AEZ","AFA","AFB","AFC","AFD","AFE","AFF","AFG","AFH","AFI","AFJ","AFK","AFL","AFM","AFN",
"AFO","AFP","AFQ","AFR","AFS","AFT","AFU","AFV","AFW","AFX","AFY","AFZ","AGA","AGB","AGC","AGD",
"AGE","AGF","AGG","AGH","AGI","AGJ","AGK","AGL","AGM","AGN","AGO","AGP","AGQ","AGR","AGS","AGT",
"AGU","AGV","AGW","AGX","AGY","AGZ","AHA","AHB","AHC","AHD","AHE","AHF","AHG","AHH","AHI","AHJ",
"AHK","AHL","AHM","AHN","AHO","AHP","AHQ","AHR","AHS","AHT","AHU","AHV","AHW","AHX","AHY","AHZ",
"AIA","AIB","AIC","AID","AIE","AIF","AIG","AIH","AII","AIJ","AIK","AIL","AIM","AIN","AIO","AIP",
"AIQ","AIR","AIS","AIT","AIU","AIV","AIW","AIX","AIY","AIZ","AJA","AJB","AJC","AJD","AJE","AJF",
"AJG","AJH","AJI","AJJ","AJK","AJL","AJM","AJN","AJO","AJP","AJQ","AJR","AJS","AJT","AJU","AJV",
"AJW","AJX","AJY","AJZ","AKA","AKB","AKC","AKD","AKE","AKF","AKG","AKH","AKI","AKJ","AKK","AKL",
"AKM","AKN","AKO","AKP","AKQ","AKR","AKS","AKT","AKU","AKV","AKW","AKX","AKY","AKZ","ALA","ALB",
"ALC","ALD","ALE","ALF","ALG","ALH","ALI","ALJ","ALK","ALL","ALM","ALN","ALO","ALP","ALQ","ALR",
"ALS","ALT","ALU","ALV","ALW","ALX","ALY","ALZ","AMA","AMB","AMC","AMD","AME","AMF","AMG","AMH",
"AMI","AMJ","AMK","AML","AMM","AMN","AMO","AMP","AMQ","AMR","AMS","AMT","AMU","AMV","AMW","AMX",
"AMY","AMZ","ANA","ANB","ANC","AND","ANE","ANF","ANG","ANH","ANI","ANJ","ANK","ANL","ANM","ANN",
"ANO","ANP","ANQ","ANR","ANS","ANT","ANU","ANV","ANW","ANX","ANY","ANZ","AOA","AOB","AOC","AOD",
"AOE","AOF","AOG","AOH","AOI","AOJ","AOK","AOL","AOM","AON","AOO","AOP","AOQ","AOR","AOS","AOT",
"AOU","AOV","AOW","AOX","AOY","AOZ","APA","APB","APC","APD","APE","APF","APG","APH","API","APJ",
"APK","APL","APM","APN","APO","APP","APQ","APR","APS","APT","APU","APV","APW","APX","APY","APZ",
"AQA","AQB","AQC","AQD","AQE","AQF","AQG","AQH","AQI","AQJ","AQK","AQL","AQM","AQN","AQO","AQP",
"AQQ","AQR","AQS","AQT","AQU","AQV","AQW","AQX","AQY","AQZ","ARA","ARB","ARC","ARD","ARE","ARF",
"ARG","ARH","ARI","ARJ","ARK","ARL","ARM","ARN","ARO","ARP","ARQ","ARR","ARS","ART","ARU","ARV",
"ARW","ARX","ARY","ARZ","ASA","ASB","ASC","ASD","ASE","ASF","ASG","ASH","ASI","ASJ","ASK","ASL",
"ASM","ASN","ASO","ASP","ASQ","ASR","ASS","AST","ASU","ASV","ASW","ASX","ASY","ASZ","ATA","ATB",
"ATC","ATD","ATE","ATF","ATG","ATH","ATI","ATJ","ATK","ATL","ATM","ATN","ATO","ATP","ATQ","ATR",
"ATS","ATT","ATU","ATV","ATW","ATX","ATY","ATZ","AUA","AUB","AUC","AUD","AUE","AUF","AUG","AUH",
"AUI","AUJ","AUK","AUL","AUM","AUN","AUO","AUP","AUQ","AUR","AUS","AUT","AUU","AUV","AUW","AUX",
"AUY","AUZ","AVA","AVB","AVC","AVD","AVE","AVF","AVG","AVH","AVI","AVJ","AVK","AVL","AVM","AVN",
"AVO","AVP","AVQ","AVR","AVS","AVT","AVU","AVV","AVW","AVX","AVY","AVZ","AWA","AWB","AWC","AWD",
"AWE","AWF","AWG","AWH","AWI","AWJ","AWK","AWL","AWM","AWN","AWO","AWP","AWQ","AWR","AWS","AWT",
"AWU","AWV","AWW","AWX","AWY","AWZ","AXA","AXB","AXC","AXD","AXE","AXF","AXG","AXH","AXI","AXJ",
"AXK","AXL","AXM","AXN","AXO","AXP","AXQ","AXR","AXS","AXT","AXU","AXV","AXW","AXX","AXY","AXZ",
"AYA","AYB","AYC","AYD","AYE","AYF","AYG","AYH","AYI","AYJ","AYK","AYL","AYM","AYN","AYO","AYP",
"AYQ","AYR","AYS","AYT","AYU","AYV","AYW","AYX","AYY","AYZ","AZA","AZB","AZC","AZD","AZE","AZF",
"AZG","AZH","AZI","AZJ","AZK","AZL","AZM","AZN","AZO","AZP","AZQ","AZR","AZS","AZT","AZU","AZV",
"AZW","AZX","AZY","AZZ","BAA","BAB","BAC","BAD","BAE","BAF","BAG","BAH","BAI","BAJ","BAK","BAL",
"BAM","BAN","BAO","BAP","BAQ","BAR","BAS","BAT","BAU","BAV","BAW","BAX","BAY","BAZ","BBA","BBB",
"BBC","BBD","BBE","BBF","BBG","BBH","BBI","BBJ","BBK","BBL","BBM","BBN","BBO","BBP","BBQ","BBR",
"BBS","BBT","BBU","BBV","BBW","BBX","BBY","BBZ","BCA","BCB","BCC","BCD","BCE","BCF","BCG","BCH",
"BCI","BCJ","BCK","BCL","BCM","BCN","BCO","BCP","BCQ","BCR","BCS","BCT","BCU","BCV","BCW","BCX",
"BCY","BCZ","BDA","BDB","BDC","BDD","BDE","BDF","BDG","BDH","BDI","BDJ","BDK","BDL","BDM","BDN",
"BDO","BDP","BDQ","BDR","BDS","BDT","BDU","BDV","BDW","BDX","BDY","BDZ","BEA","BEB","BEC","BED",
"BEE","BEF","BEG","BEH","BEI","BEJ","BEK","BEL","BEM","BEN","BEO","BEP","BEQ","BER","BES","BET",
"BEU","BEV","BEW","BEX","BEY","BEZ","BFA","BFB","BFC","BFD","BFE","BFF","BFG","BFH","BFI","BFJ",
"BFK","BFL","BFM","BFN","BFO","BFP","BFQ","BFR","BFS","BFT","BFU","BFV","BFW","BFX","BFY","BFZ",
"BGA","BGB","BGC","BGD","BGE","BGF","BGG","BGH","BGI","BGJ","BGK","BGL","BGM","BGN","BGO","BGP",
"BGQ","BGR","BGS","BGT","BGU","BGV","BGW","BGX","BGY","BGZ","BHA","BHB","BHC","BHD","BHE","BHF",
"BHG","BHH","BHI","BHJ","BHK","BHL","BHM","BHN","BHO","BHP","BHQ","BHR","BHS","BHT","BHU","BHV",
"BHW","BHX","BHY","BHZ","BIA","BIB","BIC","BID","BIE","BIF","BIG","BIH","BII","BIJ","BIK","BIL",
"BIM","BIN","BIO","BIP","BIQ","BIR","BIS","BIT","BIU","BIV","BIW","BIX","BIY","BIZ","BJA","BJB",
"BJC","BJD","BJE","BJF","BJG","BJH","BJI","BJJ","BJK","BJL","BJM","BJN","BJO","BJP","BJQ","BJR",
"BJS","BJT","BJU","BJV","BJW","BJX","BJY","BJZ","BKA","BKB","BKC","BKD","BKE","BKF","BKG","BKH",
"BKI","BKJ","BKK","BKL","BKM","BKN","BKO","BKP","BKQ","BKR","BKS","BKT","BKU","BKV","BKW","BKX",
"BKY","BKZ","BLA","BLB","BLC","BLD","BLE","BLF","BLG","BLH","BLI","BLJ","BLK","BLL","BLM","BLN",
"BLO","BLP","BLQ","BLR","BLS","BLT","BLU","BLV","BLW","BLX","BLY","BLZ","BMA","BMB","BMC","BMD",
"BME","BMF","BMG","BMH","BMI","BMJ","BMK","BML","BMM","BMN","BMO","BMP","BMQ","BMR","BMS","BMT",
"BMU","BMV","BMW","BMX","BMY","BMZ","BNA","BNB","BNC","BND","BNE","BNF","BNG","BNH","BNI","BNJ",
"BNK","BNL","BNM","BNN","BNO","BNP","BNQ","BNR","BNS","BNT","BNU","BNV","BNW","BNX","BNY","BNZ",
"BOA","BOB","BOC","BOD","BOE","BOF","BOG","BOH","BOI","BOJ","BOK","BOL","BOM","BON","BOO","BOP",
"BOQ","BOR","BOS","BOT","BOU","BOV","BOW","BOX","BOY","BOZ","BPA","BPB","BPC","BPD","BPE","BPF",
"BPG","BPH","BPI","BPJ","BPK","BPL","BPM","BPN","BPO","BPP","BPQ","BPR","BPS","BPT","BPU","BPV",
"BPW","BPX","BPY","BPZ","BQA","BQB","BQC","BQD","BQE","BQF","BQG","BQH","BQI","BQJ","BQK","BQL",
"BQM","BQN","BQO","BQP","BQQ","BQR","BQS","BQT","BQU","BQV","BQW","BQX","BQY","BQZ","BRA","BRB",
"BRC","BRD","BRE","BRF","BRG","BRH","BRI","BRJ","BRK","BRL","BRM","BRN","BRO","BRP","BRQ","BRR",
"BRS","BRT","BRU","BRV","BRW","BRX","BRY","BRZ","BSA","BSB","BSC","BSD","BSE","BSF","BSG","BSH",
"BSI","BSJ","BSK","BSL","BSM","BSN","BSO","BSP","BSQ","BSR","BSS","BST","BSU","BSV","BSW","BSX",
"BSY","BSZ","BTA","BTB","BTC","BTD","BTE","BTF","BTG","BTH","BTI","BTJ","BTK","BTL","BTM","BTN",
"BTO","BTP","BTQ","BTR","BTS","BTT","BTU","BTV","BTW","BTX","BTY","BTZ","BUA","BUB","BUC","BUD",
"BUE","BUF","BUG","BUH","BUI","BUJ","BUK","BUL","BUM","BUN","BUO","BUP","BUQ","BUR","BUS","BUT",
"BUU","BUV","BUW","BUX","BUY","BUZ","BVA","BVB","BVC","BVD","BVE","BVF","BVG","BVH","BVI","BVJ",
"BVK","BVL","BVM","BVN","BVO","BVP","BVQ","BVR","BVS","BVT","BVU","BVV","BVW","BVX","BVY","BVZ",
"BWA","BWB","BWC","BWD","BWE","BWF","BWG","BWH","BWI","BWJ","BWK","BWL","BWM","BWN","BWO","BWP",
"BWQ","BWR","BWS","BWT","BWU","BWV","BWW","BWX","BWY","BWZ","BXA","BXB","BXC","BXD","BXE","BXF",
"BXG","BXH","BXI","BXJ","BXK","BXL","BXM","BXN","BXO","BXP","BXQ","BXR","BXS","BXT","BXU","BXV",
"BXW","BXX","BXY","BXZ","BYA","BYB","BYC","BYD","BYE","BYF","BYG","BYH","BYI","BYJ","BYK","BYL",
"BYM","BYN","BYO","BYP","BYQ","BYR","BYS","BYT","BYU","BYV","BYW","BYX","BYY","BYZ","BZA","BZB",
"BZC","BZD","BZE","BZF","BZG","BZH","BZI","BZJ","BZK","BZL","BZM","BZN","BZO","BZP","BZQ","BZR",
"BZS","BZT","BZU","BZV","BZW","BZX","BZY","BZZ","CAA","CAB","CAC","CAD","CAE","CAF","CAG","CAH",
"CAI","CAJ","CAK","CAL","CAM","CAN","CAO","CAP","CAQ","CAR","CAS","CAT","CAU","CAV","CAW","CAX",
"CAY","CAZ","CBA","CBB","CBC","CBD","CBE","CBF","CBG","CBH","CBI","CBJ","CBK","CBL","CBM","CBN",
"CBO","CBP","CBQ","CBR","CBS","CBT","CBU","CBV","CBW","CBX","CBY","CBZ","CCA","CCB","CCC","CCD",
"CCE","CCF","CCG","CCH","CCI","CCJ","CCK","CCL","CCM","CCN","CCO","CCP","CCQ","CCR","CCS","CCT",
"CCU","CCV","CCW","CCX","CCY","CCZ","CDA","CDB","CDC","CDD","CDE","CDF","CDG","CDH","CDI","CDJ",
"CDK","CDL","CDM","CDN","CDO","CDP","CDQ","CDR","CDS","CDT","CDU","CDV","CDW","CDX","CDY","CDZ",
"CEA","CEB","CEC","CED","CEE","CEF","CEG","CEH","CEI","CEJ","CEK","CEL","CEM","CEN","CEO","CEP",
"CEQ","CER","CES","CET","CEU","CEV","CEW","CEX","CEY","CEZ","CFA","CFB","CFC","CFD","CFE","CFF",
"CFG","CFH","CFI","CFJ","CFK","CFL","CFM","CFN","CFO","CFP","CFQ","CFR","CFS","CFT","CFU","CFV",
"CFW","CFX","CFY","CFZ","CGA","CGB","CGC","CGD","CGE","CGF","CGG","CGH","CGI","CGJ","CGK","CGL",
"CGM","CGN","CGO","CGP","CGQ","CGR","CGS","CGT","CGU","CGV","CGW","CGX","CGY","CGZ","CHA","CHB",
"CHC","CHD","CHE","CHF","CHG","CHH","CHI","CHJ","CHK","CHL","CHM","CHN","CHO","CHP","CHQ","CHR",
"CHS","CHT","CHU","CHV","CHW","CHX","CHY","CHZ","CIA","CIB","CIC","CID","CIE","CIF","CIG","CIH",
"CII","CIJ","CIK","CIL","CIM","CIN","CIO","CIP","CIQ","CIR","CIS","CIT","CIU","CIV","CIW","CIX",
"CIY","CIZ","CJA","CJB","CJC","CJD","CJE","CJF","CJG","CJH","CJI","CJJ","CJK","CJL","CJM","CJN",
"CJO","CJP","CJQ","CJR","CJS","CJT","CJU","CJV","CJW","CJX","CJY","CJZ","CKA","CKB","CKC","CKD",
"CKE","CKF","CKG","CKH","CKI","CKJ","CKK","CKL","CKM","CKN","CKO","CKP","CKQ","CKR","CKS","CKT",
"CKU","CKV","CKW","CKX","CKY","CKZ","CLA","CLB","CLC","CLD","CLE","CLF","CLG","CLH","CLI","CLJ",
"CLK","CLL","CLM","CLN","CLO","CLP","CLQ","CLR","CLS","CLT","CLU","CLV","CLW","CLX","CLY","CLZ",
"CMA","CMB","CMC","CMD","CME","CMF","CMG","CMH","CMI","CMJ","CMK","CML","CMM","CMN","CMO","CMP",
"CMQ","CMR","CMS","CMT","CMU","CMV","CMW","CMX","CMY","CMZ","CNA","CNB","CNC","CND","CNE","CNF",
"CNG","CNH","CNI","CNJ","CNK","CNL","CNM","CNN","CNO","CNP","CNQ","CNR","CNS","CNT","CNU","CNV",
"CNW","CNX","CNY","CNZ","COA","COB","COC","COD","COE","COF","COG","COH","COI","COJ","COK","COL",
"COM","CON","COO","COP","COQ","COR","COS","COT","COU","COV","COW","COX","COY","COZ","CPA","CPB",
"CPC","CPD","CPE","CPF","CPG","CPH","CPI","CPJ","CPK","CPL","CPM","CPN","CPO","CPP","CPQ","CPR",
"CPS","CPT","CPU","CPV","CPW","CPX","CPY","CPZ","CQA","CQB","CQC","CQD","CQE","CQF","CQG","CQH",
"CQI","CQJ","CQK","CQL","CQM","CQN","CQO","CQP","CQQ","CQR","CQS","CQT","CQU","CQV","CQW","CQX",
"CQY","CQZ","CRA","CRB","CRC","CRD","CRE","CRF","CRG","CRH","CRI","CRJ","CRK","CRL","CRM","CRN",
"CRO","CRP","CRQ","CRR","CRS","CRT","CRU","CRV","CRW","CRX","CRY","CRZ","CSA","CSB","CSC","CSD",
"CSE","CSF","CSG","CSH","CSI","CSJ","CSK","CSL","CSM","CSN","CSO","CSP","CSQ","CSR","CSS","CST",
"CSU","CSV","CSW","CSX","CSY","CSZ","CTA","CTB","CTC","CTD","CTE","CTF","CTG","CTH","CTI","CTJ",
"CTK","CTL","CTM","CTN","CTO","CTP","CTQ","CTR","CTS","CTT","CTU","CTV","CTW","CTX","CTY","CTZ",
"CUA","CUB","CUC","CUD","CUE","CUF","CUG","CUH","CUI","CUJ","CUK","CUL","CUM","CUN","CUO","CUP",
"CUQ","CUR","CUS","CUT","CUU","CUV","CUW","CUX","CUY","CUZ","CVA","CVB","CVC","CVD","CVE","CVF",
"CVG","CVH","CVI","CVJ","CVK","CVL","CVM","CVN","CVO","CVP","CVQ","CVR","CVS","CVT","CVU","CVV",
"CVW","CVX","CVY","CVZ","CWA","CWB","CWC","CWD","CWE","CWF","CWG","CWH","CWI","CWJ","CWK","CWL",
"CWM","CWN","CWO","CWP","CWQ","CWR","CWS","CWT","CWU","CWV","CWW","CWX","CWY","CWZ","CXA","CXB",
"CXC","CXD","CXE","CXF","CXG","CXH","CXI","CXJ","CXK","CXL","CXM","CXN","CXO","CXP","CXQ","CXR",
"CXS","CXT","CXU","CXV","CXW","CXX","CXY","CXZ","CYA","CYB","CYC","CYD","CYE","CYF","CYG","CYH",
"CYI","CYJ","CYK","CYL","CYM","CYN","CYO","CYP","CYQ","CYR","CYS","CYT","CYU","CYV","CYW","CYX",
"CYY","CYZ","CZA","CZB","CZC","CZD","CZE","CZF","CZG","CZH","CZI","CZJ","CZK","CZL","CZM","CZN",
"CZO","CZP","CZQ","CZR","CZS","CZT","CZU","CZV","CZW","CZX","CZY","CZZ","DAA","DAB","DAC","DAD",
"DAE","DAF","DAG","DAH","DAI","DAJ","DAK","DAL","DAM","DAN","DAO","DAP","DAQ","DAR","DAS","DAT",
"DAU","DAV","DAW","DAX","DAY","DAZ","DBA","DBB","DBC","DBD","DBE","DBF","DBG","DBH","DBI","DBJ",
"DBK","DBL","DBM","DBN","DBO","DBP","DBQ","DBR","DBS","DBT","DBU","DBV","DBW","DBX","DBY","DBZ",
"DCA","DCB","DCC","DCD","DCE","DCF","DCG","DCH","DCI","DCJ","DCK","DCL","DCM","DCN","DCO","DCP",
"DCQ","DCR","DCS","DCT","DCU","DCV","DCW","DCX","DCY","DCZ","DDA","DDB","DDC","DDD","DDE","DDF",
"DDG","DDH","DDI","DDJ","DDK","DDL","DDM","DDN","DDO","DDP","DDQ","DDR","DDS","DDT","DDU","DDV",
"DDW","DDX","DDY","DDZ","DEA","DEB","DEC","DED","DEE","DEF","DEG","DEH","DEI","DEJ","DEK","DEL",
"DEM","DEN","DEO","DEP","DEQ","DER","DES","DET","DEU","DEV","DEW","DEX","DEY","DEZ","DFA","DFB",
"DFC","DFD","DFE","DFF","DFG","DFH","DFI","DFJ","DFK","DFL","DFM","DFN","DFO","DFP","DFQ","DFR",
"DFS","DFT","DFU","DFV","DFW","DFX","DFY","DFZ","DGA","DGB","DGC","DGD","DGE","DGF","DGG","DGH",
"DGI","DGJ","DGK","DGL","DGM","DGN","DGO","DGP","DGQ","DGR","DGS","DGT","DGU","DGV","DGW","DGX",
"DGY","DGZ","DHA","DHB","DHC","DHD","DHE","DHF","DHG","DHH","DHI","DHJ","DHK","DHL","DHM","DHN",
"DHO","DHP","DHQ","DHR","DHS","DHT","DHU","DHV","DHW","DHX","DHY","DHZ","DIA","DIB","DIC","DID",
"DIE","DIF","DIG","DIH","DII","DIJ","DIK","DIL","DIM","DIN","DIO","DIP","DIQ","DIR","DIS","DIT",
"DIU","DIV","DIW","DIX","DIY","DIZ","DJA","DJB","DJC","DJD","DJE","DJF","DJG","DJH","DJI","DJJ",
"DJK","DJL","DJM","DJN","DJO","DJP","DJQ","DJR","DJS","DJT","DJU","DJV","DJW","DJX","DJY","DJZ",
"DKA","DKB","DKC","DKD","DKE","DKF","DKG","DKH","DKI","DKJ","DKK","DKL","DKM","DKN","DKO","DKP",
"DKQ","DKR","DKS","DKT","DKU","DKV","DKW","DKX","DKY","DKZ","DLA","DLB","DLC","DLD","DLE","DLF",
"DLG","DLH","DLI","DLJ","DLK","DLL","DLM","DLN","DLO","DLP","DLQ","DLR","DLS","DLT","DLU","DLV",
"DLW","DLX","DLY","DLZ","DMA","DMB","DMC","DMD","DME","DMF","DMG","DMH","DMI","DMJ","DMK","DML",
"DMM","DMN","DMO","DMP","DMQ","DMR","DMS","DMT","DMU","DMV","DMW","DMX","DMY","DMZ","DNA","DNB",
"DNC","DND","DNE","DNF","DNG","DNH","DNI","DNJ","DNK","DNL","DNM","DNN","DNO","DNP","DNQ","DNR",
"DNS","DNT","DNU","DNV","DNW","DNX","DNY","DNZ","DOA","DOB","DOC","DOD","DOE","DOF","DOG","DOH",
"DOI","DOJ","DOK","DOL","DOM","DON","DOO","DOP","DOQ","DOR","DOS","DOT","DOU","DOV","DOW","DOX",
"DOY","DOZ","DPA","DPB","DPC","DPD","DPE","DPF","DPG","DPH","DPI","DPJ","DPK","DPL","DPM","DPN",
"DPO","DPP","DPQ","DPR","DPS","DPT","DPU","DPV","DPW","DPX","DPY","DPZ","DQA","DQB","DQC","DQD",
"DQE","DQF","DQG","DQH","DQI","DQJ","DQK","DQL","DQM","DQN","DQO","DQP","DQQ","DQR","DQS","DQT",
"DQU","DQV","DQW","DQX","DQY","DQZ","DRA","DRB","DRC","DRD","DRE","DRF","DRG","DRH","DRI","DRJ",
"DRK","DRL","DRM","DRN","DRO","DRP","DRQ","DRR","DRS","DRT","DRU","DRV","DRW","DRX","DRY","DRZ",
"DSA","DSB","DSC","DSD","DSE","DSF","DSG","DSH","DSI","DSJ","DSK","DSL","DSM","DSN","DSO","DSP",
"DSQ","DSR","DSS","DST","DSU","DSV","DSW","DSX","DSY","DSZ","DTA","DTB","DTC","DTD","DTE","DTF",
"DTG","DTH","DTI","DTJ","DTK","DTL","DTM","DTN","DTO","DTP","DTQ","DTR","DTS","DTT","DTU","DTV",
"DTW","DTX","DTY","DTZ","DUA","DUB","DUC","DUD","DUE","DUF","DUG","DUH","DUI","DUJ","DUK","DUL",
"DUM","DUN","DUO","DUP","DUQ","DUR","DUS","DUT","DUU","DUV","DUW","DUX","DUY","DUZ","DVA","DVB",
"DVC","DVD","DVE","DVF","DVG","DVH","DVI","DVJ","DVK","DVL","DVM","DVN","DVO","DVP","DVQ","DVR",
"DVS","DVT","DVU","DVV","DVW","DVX","DVY","DVZ","DWA","DWB","DWC","DWD","DWE","DWF","DWG","DWH",
"DWI","DWJ","DWK","DWL","DWM","DWN","DWO","DWP","DWQ","DWR","DWS","DWT","DWU","DWV","DWW","DWX",
"DWY","DWZ","DXA","DXB","DXC","DXD","DXE","DXF","DXG","DXH","DXI","DXJ","DXK","DXL","DXM","DXN",
"DXO","DXP","DXQ","DXR","DXS","DXT","DXU","DXV","DXW","DXX","DXY","DXZ","DYA","DYB","DYC","DYD",
"DYE","DYF","DYG","DYH","DYI","DYJ","DYK","DYL","DYM","DYN","DYO","DYP","DYQ","DYR","DYS","DYT",
"DYU","DYV","DYW","DYX","DYY","DYZ","DZA","DZB","DZC","DZD","DZE","DZF","DZG","DZH","DZI","DZJ",
"DZK","DZL","DZM","DZN","DZO","DZP","DZQ","DZR","DZS","DZT","DZU","DZV","DZW","DZX","DZY","DZZ",
/*
    E-started's intentionally omitted
*/
"FAA","FAB","FAC","FAD","FAE","FAF","FAG","FAH","FAI","FAJ","FAK","FAL",
"FAM","FAN","FAO","FAP","FAQ","FAR","FAS","FAT","FAU","FAV","FAW","FAX","FAY","FAZ","FBA","FBB",
"FBC","FBD","FBE","FBF","FBG","FBH","FBI","FBJ","FBK","FBL","FBM","FBN","FBO","FBP","FBQ","FBR",
"FBS","FBT","FBU","FBV","FBW","FBX","FBY","FBZ","FCA","FCB","FCC","FCD","FCE","FCF","FCG","FCH",
"FCI","FCJ","FCK","FCL","FCM","FCN","FCO","FCP","FCQ","FCR","FCS","FCT","FCU","FCV","FCW","FCX",
"FCY","FCZ","FDA","FDB","FDC","FDD","FDE","FDF","FDG","FDH","FDI","FDJ","FDK","FDL","FDM","FDN",
"FDO","FDP","FDQ","FDR","FDS","FDT","FDU","FDV","FDW","FDX","FDY","FDZ","FEA","FEB","FEC","FED",
"FEE","FEF","FEG","FEH","FEI","FEJ","FEK","FEL","FEM","FEN","FEO","FEP","FEQ","FER","FES","FET",
"FEU","FEV","FEW","FEX","FEY","FEZ","FFA","FFB","FFC","FFD","FFE","FFF","FFG","FFH","FFI","FFJ",
"FFK","FFL","FFM","FFN","FFO","FFP","FFQ","FFR","FFS","FFT","FFU","FFV","FFW","FFX","FFY","FFZ",
"FGA","FGB","FGC","FGD","FGE","FGF","FGG","FGH","FGI","FGJ","FGK","FGL","FGM","FGN","FGO","FGP",
"FGQ","FGR","FGS","FGT","FGU","FGV","FGW","FGX","FGY","FGZ","FHA","FHB","FHC","FHD","FHE","FHF",
"FHG","FHH","FHI","FHJ","FHK","FHL","FHM","FHN","FHO","FHP","FHQ","FHR","FHS","FHT","FHU","FHV",
"FHW","FHX","FHY","FHZ","FIA","FIB","FIC","FID","FIE","FIF","FIG","FIH","FII","FIJ","FIK","FIL",
"FIM","FIN","FIO","FIP","FIQ","FIR","FIS","FIT","FIU","FIV","FIW","FIX","FIY","FIZ","FJA","FJB",
"FJC","FJD","FJE","FJF","FJG","FJH","FJI","FJJ","FJK","FJL","FJM","FJN","FJO","FJP","FJQ","FJR",
"FJS","FJT","FJU","FJV","FJW","FJX","FJY","FJZ","FKA","FKB","FKC","FKD","FKE","FKF","FKG","FKH",
"FKI","FKJ","FKK","FKL","FKM","FKN","FKO","FKP","FKQ","FKR","FKS","FKT","FKU","FKV","FKW","FKX",
"FKY","FKZ","FLA","FLB","FLC","FLD","FLE","FLF","FLG","FLH","FLI","FLJ","FLK","FLL","FLM","FLN",
"FLO","FLP","FLQ","FLR","FLS","FLT","FLU","FLV","FLW","FLX","FLY","FLZ","FMA","FMB","FMC","FMD",
"FME","FMF","FMG","FMH","FMI","FMJ","FMK","FML","FMM","FMN","FMO","FMP","FMQ","FMR","FMS","FMT",
"FMU","FMV","FMW","FMX","FMY","FMZ","FNA","FNB","FNC","FND","FNE","FNF","FNG","FNH","FNI","FNJ",
"FNK","FNL","FNM","FNN","FNO","FNP","FNQ","FNR","FNS","FNT","FNU","FNV","FNW","FNX","FNY","FNZ",
"FOA","FOB","FOC","FOD","FOE","FOF","FOG","FOH","FOI","FOJ","FOK","FOL","FOM","FON","FOO","FOP",
"FOQ","FOR","FOS","FOT","FOU","FOV","FOW","FOX","FOY","FOZ","FPA","FPB","FPC","FPD","FPE","FPF",
"FPG","FPH","FPI","FPJ","FPK","FPL","FPM","FPN","FPO","FPP","FPQ","FPR","FPS","FPT","FPU","FPV",
"FPW","FPX","FPY","FPZ","FQA","FQB","FQC","FQD","FQE","FQF","FQG","FQH","FQI","FQJ","FQK","FQL",
"FQM","FQN","FQO","FQP","FQQ","FQR","FQS","FQT","FQU","FQV","FQW","FQX","FQY","FQZ","FRA","FRB",
"FRC","FRD","FRE","FRF","FRG","FRH","FRI","FRJ","FRK","FRL","FRM","FRN","FRO","FRP","FRQ","FRR",
"FRS","FRT","FRU","FRV","FRW","FRX","FRY","FRZ","FSA","FSB","FSC","FSD","FSE","FSF","FSG","FSH",
"FSI","FSJ","FSK","FSL","FSM","FSN","FSO","FSP","FSQ","FSR","FSS","FST","FSU","FSV","FSW","FSX",
"FSY","FSZ","FTA","FTB","FTC","FTD","FTE","FTF","FTG","FTH","FTI","FTJ","FTK","FTL","FTM","FTN",
"FTO","FTP","FTQ","FTR","FTS","FTT","FTU","FTV","FTW","FTX","FTY","FTZ","FUA","FUB","FUC","FUD",
"FUE","FUF","FUG","FUH","FUI","FUJ","FUK","FUL","FUM","FUN","FUO","FUP","FUQ","FUR","FUS","FUT",
"FUU","FUV","FUW","FUX","FUY","FUZ","FVA","FVB","FVC","FVD","FVE","FVF","FVG","FVH","FVI","FVJ",
"FVK","FVL","FVM","FVN","FVO","FVP","FVQ","FVR","FVS","FVT","FVU","FVV","FVW","FVX","FVY","FVZ",
"FWA","FWB","FWC","FWD","FWE","FWF","FWG","FWH","FWI","FWJ","FWK","FWL","FWM","FWN","FWO","FWP",
"FWQ","FWR","FWS","FWT","FWU","FWV","FWW","FWX","FWY","FWZ","FXA","FXB","FXC","FXD","FXE","FXF",
"FXG","FXH","FXI","FXJ","FXK","FXL","FXM","FXN","FXO","FXP","FXQ","FXR","FXS","FXT","FXU","FXV",
"FXW","FXX","FXY","FXZ","FYA","FYB","FYC","FYD","FYE","FYF","FYG","FYH","FYI","FYJ","FYK","FYL",
"FYM","FYN","FYO","FYP","FYQ","FYR","FYS","FYT","FYU","FYV","FYW","FYX","FYY","FYZ","FZA","FZB",
"FZC","FZD","FZE","FZF","FZG","FZH","FZI","FZJ","FZK","FZL","FZM","FZN","FZO","FZP","FZQ","FZR",
"FZS","FZT","FZU","FZV","FZW","FZX","FZY","FZZ","GAA","GAB","GAC","GAD","GAE","GAF","GAG","GAH",
"GAI","GAJ","GAK","GAL","GAM","GAN","GAO","GAP","GAQ","GAR","GAS","GAT","GAU","GAV","GAW","GAX",
"GAY","GAZ","GBA","GBB","GBC","GBD","GBE","GBF","GBG","GBH","GBI","GBJ","GBK","GBL","GBM","GBN",
"GBO","GBP","GBQ","GBR","GBS","GBT","GBU","GBV","GBW","GBX","GBY","GBZ","GCA","GCB","GCC","GCD",
"GCE","GCF","GCG","GCH","GCI","GCJ","GCK","GCL","GCM","GCN","GCO","GCP","GCQ","GCR","GCS","GCT",
"GCU","GCV","GCW","GCX","GCY","GCZ","GDA","GDB","GDC","GDD","GDE","GDF","GDG","GDH","GDI","GDJ",
"GDK","GDL","GDM","GDN","GDO","GDP","GDQ","GDR","GDS","GDT","GDU","GDV","GDW","GDX","GDY","GDZ",
"GEA","GEB","GEC","GED","GEE","GEF","GEG","GEH","GEI","GEJ","GEK","GEL","GEM","GEN","GEO","GEP",
"GEQ","GER","GES","GET","GEU","GEV","GEW","GEX","GEY","GEZ","GFA","GFB","GFC","GFD","GFE","GFF",
"GFG","GFH","GFI","GFJ","GFK","GFL","GFM","GFN","GFO","GFP","GFQ","GFR","GFS","GFT","GFU","GFV",
"GFW","GFX","GFY","GFZ","GGA","GGB","GGC","GGD","GGE","GGF","GGG","GGH","GGI","GGJ","GGK","GGL",
"GGM","GGN","GGO","GGP","GGQ","GGR","GGS","GGT","GGU","GGV","GGW","GGX","GGY","GGZ","GHA","GHB",
"GHC","GHD","GHE","GHF","GHG","GHH","GHI","GHJ","GHK","GHL","GHM","GHN","GHO","GHP","GHQ","GHR",
"GHS","GHT","GHU","GHV","GHW","GHX","GHY","GHZ","GIA","GIB","GIC","GID","GIE","GIF","GIG","GIH",
"GII","GIJ","GIK","GIL","GIM","GIN","GIO","GIP","GIQ","GIR","GIS","GIT","GIU","GIV","GIW","GIX",
"GIY","GIZ","GJA","GJB","GJC","GJD","GJE","GJF","GJG","GJH","GJI","GJJ","GJK","GJL","GJM","GJN",
"GJO","GJP","GJQ","GJR","GJS","GJT","GJU","GJV","GJW","GJX","GJY","GJZ","GKA","GKB","GKC","GKD",
"GKE","GKF","GKG","GKH","GKI","GKJ","GKK","GKL","GKM","GKN","GKO","GKP","GKQ","GKR","GKS","GKT",
"GKU","GKV","GKW","GKX","GKY","GKZ","GLA","GLB","GLC","GLD","GLE","GLF","GLG","GLH","GLI","GLJ",
"GLK","GLL","GLM","GLN","GLO","GLP","GLQ","GLR","GLS","GLT","GLU","GLV","GLW","GLX","GLY","GLZ",
"GMA","GMB","GMC","GMD","GME","GMF","GMG","GMH","GMI","GMJ","GMK","GML","GMM","GMN","GMO","GMP",
"GMQ","GMR","GMS","GMT","GMU","GMV","GMW","GMX","GMY","GMZ","GNA","GNB","GNC","GND","GNE","GNF",
"GNG","GNH","GNI","GNJ","GNK","GNL","GNM","GNN","GNO","GNP","GNQ","GNR","GNS","GNT","GNU","GNV",
"GNW","GNX","GNY","GNZ","GOA","GOB","GOC","GOD","GOE","GOF","GOG","GOH","GOI","GOJ","GOK","GOL",
"GOM","GON","GOO","GOP","GOQ","GOR","GOS","GOT","GOU","GOV","GOW","GOX","GOY","GOZ","GPA","GPB",
"GPC","GPD","GPE","GPF","GPG","GPH","GPI","GPJ","GPK","GPL","GPM","GPN","GPO","GPP","GPQ","GPR",
"GPS","GPT","GPU","GPV","GPW","GPX","GPY","GPZ","GQA","GQB","GQC","GQD","GQE","GQF","GQG","GQH",
"GQI","GQJ","GQK","GQL","GQM","GQN","GQO","GQP","GQQ","GQR","GQS","GQT","GQU","GQV","GQW","GQX",
"GQY","GQZ","GRA","GRB","GRC","GRD","GRE","GRF","GRG","GRH","GRI","GRJ","GRK","GRL","GRM","GRN",
"GRO","GRP","GRQ","GRR","GRS","GRT","GRU","GRV","GRW","GRX","GRY","GRZ","GSA","GSB","GSC","GSD",
"GSE","GSF","GSG","GSH","GSI","GSJ","GSK","GSL","GSM","GSN","GSO","GSP","GSQ","GSR","GSS","GST",
"GSU","GSV","GSW","GSX","GSY","GSZ","GTA","GTB","GTC","GTD","GTE","GTF","GTG","GTH","GTI","GTJ",
"GTK","GTL","GTM","GTN","GTO","GTP","GTQ","GTR","GTS","GTT","GTU","GTV","GTW","GTX","GTY","GTZ",
"GUA","GUB","GUC","GUD","GUE","GUF","GUG","GUH","GUI","GUJ","GUK","GUL","GUM","GUN","GUO","GUP",
"GUQ","GUR","GUS","GUT","GUU","GUV","GUW","GUX","GUY","GUZ","GVA","GVB","GVC","GVD","GVE","GVF",
"GVG","GVH","GVI","GVJ","GVK","GVL","GVM","GVN","GVO","GVP","GVQ","GVR","GVS","GVT","GVU","GVV",
"GVW","GVX","GVY","GVZ","GWA","GWB","GWC","GWD","GWE","GWF","GWG","GWH","GWI","GWJ","GWK","GWL",
"GWM","GWN","GWO","GWP","GWQ","GWR","GWS","GWT","GWU","GWV","GWW","GWX","GWY","GWZ","GXA","GXB",
"GXC","GXD","GXE","GXF","GXG","GXH","GXI","GXJ","GXK","GXL","GXM","GXN","GXO","GXP","GXQ","GXR",
"GXS","GXT","GXU","GXV","GXW","GXX","GXY","GXZ","GYA","GYB","GYC","GYD","GYE","GYF","GYG","GYH",
"GYI","GYJ","GYK","GYL","GYM","GYN","GYO","GYP","GYQ","GYR","GYS","GYT","GYU","GYV","GYW","GYX",
"GYY","GYZ","GZA","GZB","GZC","GZD","GZE","GZF","GZG","GZH","GZI","GZJ","GZK","GZL","GZM","GZN",
"GZO","GZP","GZQ","GZR","GZS","GZT","GZU","GZV","GZW","GZX","GZY","GZZ","HAA","HAB","HAC","HAD",
"HAE","HAF","HAG","HAH","HAI","HAJ","HAK","HAL","HAM","HAN","HAO","HAP","HAQ","HAR","HAS","HAT",
"HAU","HAV","HAW","HAX","HAY","HAZ","HBA","HBB","HBC","HBD","HBE","HBF","HBG","HBH","HBI","HBJ",
"HBK","HBL","HBM","HBN","HBO","HBP","HBQ","HBR","HBS","HBT","HBU","HBV","HBW","HBX","HBY","HBZ",
"HCA","HCB","HCC","HCD","HCE","HCF","HCG","HCH","HCI","HCJ","HCK","HCL","HCM","HCN","HCO","HCP",
"HCQ","HCR","HCS","HCT","HCU","HCV","HCW","HCX","HCY","HCZ","HDA","HDB","HDC","HDD","HDE","HDF",
"HDG","HDH","HDI","HDJ","HDK","HDL","HDM","HDN","HDO","HDP","HDQ","HDR","HDS","HDT","HDU","HDV",
"HDW","HDX","HDY","HDZ","HEA","HEB","HEC","HED","HEE","HEF","HEG","HEH","HEI","HEJ","HEK","HEL",
"HEM","HEN","HEO","HEP","HEQ","HER","HES","HET","HEU","HEV","HEW","HEX","HEY","HEZ","HFA","HFB",
"HFC","HFD","HFE","HFF","HFG","HFH","HFI","HFJ","HFK","HFL","HFM","HFN","HFO","HFP","HFQ","HFR",
"HFS","HFT","HFU","HFV","HFW","HFX","HFY","HFZ","HGA","HGB","HGC","HGD","HGE","HGF","HGG","HGH",
"HGI","HGJ","HGK","HGL","HGM","HGN","HGO","HGP","HGQ","HGR","HGS","HGT","HGU","HGV","HGW","HGX",
"HGY","HGZ","HHA","HHB","HHC","HHD","HHE","HHF","HHG","HHH","HHI","HHJ","HHK","HHL","HHM","HHN",
"HHO","HHP","HHQ","HHR","HHS","HHT","HHU","HHV","HHW","HHX","HHY","HHZ","HIA","HIB","HIC","HID",
"HIE","HIF","HIG","HIH","HII","HIJ","HIK","HIL","HIM","HIN","HIO","HIP","HIQ","HIR","HIS","HIT",
"HIU","HIV","HIW","HIX","HIY","HIZ","HJA","HJB","HJC","HJD","HJE","HJF","HJG","HJH","HJI","HJJ",
"HJK","HJL","HJM","HJN","HJO","HJP","HJQ","HJR","HJS","HJT","HJU","HJV","HJW","HJX","HJY","HJZ",
"HKA","HKB","HKC","HKD","HKE","HKF","HKG","HKH","HKI","HKJ","HKK","HKL","HKM","HKN","HKO","HKP",
"HKQ","HKR","HKS","HKT","HKU","HKV","HKW","HKX","HKY","HKZ","HLA","HLB","HLC","HLD","HLE","HLF",
"HLG","HLH","HLI","HLJ","HLK","HLL","HLM","HLN","HLO","HLP","HLQ","HLR","HLS","HLT","HLU","HLV",
"HLW","HLX","HLY","HLZ","HMA","HMB","HMC","HMD","HME","HMF","HMG","HMH","HMI","HMJ","HMK","HML",
"HMM","HMN","HMO","HMP","HMQ","HMR","HMS","HMT","HMU","HMV","HMW","HMX","HMY","HMZ","HNA","HNB",
"HNC","HND","HNE","HNF","HNG","HNH","HNI","HNJ","HNK","HNL","HNM","HNN","HNO","HNP","HNQ","HNR",
"HNS","HNT","HNU","HNV","HNW","HNX","HNY","HNZ","HOA","HOB","HOC","HOD","HOE","HOF","HOG","HOH",
"HOI","HOJ","HOK","HOL","HOM","HON","HOO","HOP","HOQ","HOR","HOS","HOT","HOU","HOV","HOW","HOX",
"HOY","HOZ","HPA","HPB","HPC","HPD","HPE","HPF","HPG","HPH","HPI","HPJ","HPK","HPL","HPM","HPN",
"HPO","HPP","HPQ","HPR","HPS","HPT","HPU","HPV","HPW","HPX","HPY","HPZ","HQA","HQB","HQC","HQD",
"HQE","HQF","HQG","HQH","HQI","HQJ","HQK","HQL","HQM","HQN","HQO","HQP","HQQ","HQR","HQS","HQT",
"HQU","HQV","HQW","HQX","HQY","HQZ","HRA","HRB","HRC","HRD","HRE","HRF","HRG","HRH","HRI","HRJ",
"HRK","HRL","HRM","HRN","HRO","HRP","HRQ","HRR","HRS","HRT","HRU","HRV","HRW","HRX","HRY","HRZ",
"HSA","HSB","HSC","HSD","HSE","HSF","HSG","HSH","HSI","HSJ","HSK","HSL","HSM","HSN","HSO","HSP",
"HSQ","HSR","HSS","HST","HSU","HSV","HSW","HSX","HSY","HSZ","HTA","HTB","HTC","HTD","HTE","HTF",
"HTG","HTH","HTI","HTJ","HTK","HTL","HTM","HTN","HTO","HTP","HTQ","HTR","HTS","HTT","HTU","HTV",
"HTW","HTX","HTY","HTZ","HUA","HUB","HUC","HUD","HUE","HUF","HUG","HUH","HUI","HUJ","HUK","HUL",
"HUM","HUN","HUO","HUP","HUQ","HUR","HUS","HUT","HUU","HUV","HUW","HUX","HUY","HUZ","HVA","HVB",
"HVC","HVD","HVE","HVF","HVG","HVH","HVI","HVJ","HVK","HVL","HVM","HVN","HVO","HVP","HVQ","HVR",
"HVS","HVT","HVU","HVV","HVW","HVX","HVY","HVZ","HWA","HWB","HWC","HWD","HWE","HWF","HWG","HWH",
"HWI","HWJ","HWK","HWL","HWM","HWN","HWO","HWP","HWQ","HWR","HWS","HWT","HWU","HWV","HWW","HWX",
"HWY","HWZ","HXA","HXB","HXC","HXD","HXE","HXF","HXG","HXH","HXI","HXJ","HXK","HXL","HXM","HXN",
"HXO","HXP","HXQ","HXR","HXS","HXT","HXU","HXV","HXW","HXX","HXY","HXZ","HYA","HYB","HYC","HYD",
"HYE","HYF","HYG","HYH","HYI","HYJ","HYK","HYL","HYM","HYN","HYO","HYP","HYQ","HYR","HYS","HYT",
"HYU","HYV","HYW","HYX","HYY","HYZ","HZA","HZB","HZC","HZD","HZE","HZF","HZG","HZH","HZI","HZJ",
"HZK","HZL","HZM","HZN","HZO","HZP","HZQ","HZR","HZS","HZT","HZU","HZV","HZW","HZX","HZY","HZZ",
"IAA","IAB","IAC","IAD","IAE","IAF","IAG","IAH","IAI","IAJ","IAK","IAL","IAM","IAN","IAO","IAP",
"IAQ","IAR","IAS","IAT","IAU","IAV","IAW","IAX","IAY","IAZ","IBA","IBB","IBC","IBD","IBE","IBF",
"IBG","IBH","IBI","IBJ","IBK","IBL","IBM","IBN","IBO","IBP","IBQ","IBR","IBS","IBT","IBU","IBV",
"IBW","IBX","IBY","IBZ","ICA","ICB","ICC","ICD","ICE","ICF","ICG","ICH","ICI","ICJ","ICK","ICL",
"ICM","ICN","ICO","ICP","ICQ","ICR","ICS","ICT","ICU","ICV","ICW","ICX","ICY","ICZ","IDA","IDB",
"IDC","IDD","IDE","IDF","IDG","IDH","IDI","IDJ","IDK","IDL","IDM","IDN","IDO","IDP","IDQ","IDR",
"IDS","IDT","IDU","IDV","IDW","IDX","IDY","IDZ","IEA","IEB","IEC","IED","IEE","IEF","IEG","IEH",
"IEI","IEJ","IEK","IEL","IEM","IEN","IEO","IEP","IEQ","IER","IES","IET","IEU","IEV","IEW","IEX",
"IEY","IEZ","IFA","IFB","IFC","IFD","IFE","IFF","IFG","IFH","IFI","IFJ","IFK","IFL","IFM","IFN",
"IFO","IFP","IFQ","IFR","IFS","IFT","IFU","IFV","IFW","IFX","IFY","IFZ","IGA","IGB","IGC","IGD",
"IGE","IGF","IGG","IGH","IGI","IGJ","IGK","IGL","IGM","IGN","IGO","IGP","IGQ","IGR","IGS","IGT",
"IGU","IGV","IGW","IGX","IGY","IGZ","IHA","IHB","IHC","IHD","IHE","IHF","IHG","IHH","IHI","IHJ",
"IHK","IHL","IHM","IHN","IHO","IHP","IHQ","IHR","IHS","IHT","IHU","IHV","IHW","IHX","IHY","IHZ",
"IIA","IIB","IIC","IID","IIE","IIF","IIG","IIH","III","IIJ","IIK","IIL","IIM","IIN","IIO","IIP",
"IIQ","IIR","IIS","IIT","IIU","IIV","IIW","IIX","IIY","IIZ","IJA","IJB","IJC","IJD","IJE","IJF",
"IJG","IJH","IJI","IJJ","IJK","IJL","IJM","IJN","IJO","IJP","IJQ","IJR","IJS","IJT","IJU","IJV",
"IJW","IJX","IJY","IJZ","IKA","IKB","IKC","IKD","IKE","IKF","IKG","IKH","IKI","IKJ","IKK","IKL",
"IKM","IKN","IKO","IKP","IKQ","IKR","IKS","IKT","IKU","IKV","IKW","IKX","IKY","IKZ","ILA","ILB",
"ILC","ILD","ILE","ILF","ILG","ILH","ILI","ILJ","ILK","ILL","ILM","ILN","ILO","ILP","ILQ","ILR",
"ILS","ILT","ILU","ILV","ILW","ILX","ILY","ILZ","IMA","IMB","IMC","IMD","IME","IMF","IMG","IMH",
"IMI","IMJ","IMK","IML","IMM","IMN","IMO","IMP","IMQ","IMR","IMS","IMT","IMU","IMV","IMW","IMX",
"IMY","IMZ","INA","INB","INC","IND","INE","INF","ING","INH","INI","INJ","INK","INL","INM","INN",
"INO","INP","INQ","INR","INS","INT","INU","INV","INW","INX","INY","INZ","IOA","IOB","IOC","IOD",
"IOE","IOF","IOG","IOH","IOI","IOJ","IOK","IOL","IOM","ION","IOO","IOP","IOQ","IOR","IOS","IOT",
"IOU","IOV","IOW","IOX","IOY","IOZ","IPA","IPB","IPC","IPD","IPE","IPF","IPG","IPH","IPI","IPJ",
"IPK","IPL","IPM","IPN","IPO","IPP","IPQ","IPR","IPS","IPT","IPU","IPV","IPW","IPX","IPY","IPZ",
"IQA","IQB","IQC","IQD","IQE","IQF","IQG","IQH","IQI","IQJ","IQK","IQL","IQM","IQN","IQO","IQP",
"IQQ","IQR","IQS","IQT","IQU","IQV","IQW","IQX","IQY","IQZ","IRA","IRB","IRC","IRD","IRE","IRF",
"IRG","IRH","IRI","IRJ","IRK","IRL","IRM","IRN","IRO","IRP","IRQ","IRR","IRS","IRT","IRU","IRV",
"IRW","IRX","IRY","IRZ","ISA","ISB","ISC","ISD","ISE","ISF","ISG","ISH","ISI","ISJ","ISK","ISL",
"ISM","ISN","ISO","ISP","ISQ","ISR","ISS","IST","ISU","ISV","ISW","ISX","ISY","ISZ","ITA","ITB",
"ITC","ITD","ITE","ITF","ITG","ITH","ITI","ITJ","ITK","ITL","ITM","ITN","ITO","ITP","ITQ","ITR",
"ITS","ITT","ITU","ITV","ITW","ITX","ITY","ITZ","IUA","IUB","IUC","IUD","IUE","IUF","IUG","IUH",
"IUI","IUJ","IUK","IUL","IUM","IUN","IUO","IUP","IUQ","IUR","IUS","IUT","IUU","IUV","IUW","IUX",
"IUY","IUZ","IVA","IVB","IVC","IVD","IVE","IVF","IVG","IVH","IVI","IVJ","IVK","IVL","IVM","IVN",
"IVO","IVP","IVQ","IVR","IVS","IVT","IVU","IVV","IVW","IVX","IVY","IVZ","IWA","IWB","IWC","IWD",
"IWE","IWF","IWG","IWH","IWI","IWJ","IWK","IWL","IWM","IWN","IWO","IWP","IWQ","IWR","IWS","IWT",
"IWU","IWV","IWW","IWX","IWY","IWZ","IXA","IXB","IXC","IXD","IXE","IXF","IXG","IXH","IXI","IXJ",
"IXK","IXL","IXM","IXN","IXO","IXP","IXQ","IXR","IXS","IXT","IXU","IXV","IXW","IXX","IXY","IXZ",
"IYA","IYB","IYC","IYD","IYE","IYF","IYG","IYH","IYI","IYJ","IYK","IYL","IYM","IYN","IYO","IYP",
"IYQ","IYR","IYS","IYT","IYU","IYV","IYW","IYX","IYY","IYZ","IZA","IZB","IZC","IZD","IZE","IZF",
"IZG","IZH","IZI","IZJ","IZK","IZL","IZM","IZN","IZO","IZP","IZQ","IZR","IZS","IZT","IZU","IZV",
"IZW","IZX","IZY","IZZ","JAA","JAB","JAC","JAD","JAE","JAF","JAG","JAH","JAI","JAJ","JAK","JAL",
"JAM","JAN","JAO","JAP","JAQ","JAR","JAS","JAT","JAU","JAV","JAW","JAX","JAY","JAZ","JBA","JBB",
"JBC","JBD","JBE","JBF","JBG","JBH","JBI","JBJ","JBK","JBL","JBM","JBN","JBO","JBP","JBQ","JBR",
"JBS","JBT","JBU","JBV","JBW","JBX","JBY","JBZ","JCA","JCB","JCC","JCD","JCE","JCF","JCG","JCH",
"JCI","JCJ","JCK","JCL","JCM","JCN","JCO","JCP","JCQ","JCR","JCS","JCT","JCU","JCV","JCW","JCX",
"JCY","JCZ","JDA","JDB","JDC","JDD","JDE","JDF","JDG","JDH","JDI","JDJ","JDK","JDL","JDM","JDN",
"JDO","JDP","JDQ","JDR","JDS","JDT","JDU","JDV","JDW","JDX","JDY","JDZ","JEA","JEB","JEC","JED",
"JEE","JEF","JEG","JEH","JEI","JEJ","JEK","JEL","JEM","JEN","JEO","JEP","JEQ","JER","JES","JET",
"JEU","JEV","JEW","JEX","JEY","JEZ","JFA","JFB","JFC","JFD","JFE","JFF","JFG","JFH","JFI","JFJ",
"JFK","JFL","JFM","JFN","JFO","JFP","JFQ","JFR","JFS","JFT","JFU","JFV","JFW","JFX","JFY","JFZ",
"JGA","JGB","JGC","JGD","JGE","JGF","JGG","JGH","JGI","JGJ","JGK","JGL","JGM","JGN","JGO","JGP",
"JGQ","JGR","JGS","JGT","JGU","JGV","JGW","JGX","JGY","JGZ","JHA","JHB","JHC","JHD","JHE","JHF",
"JHG","JHH","JHI","JHJ","JHK","JHL","JHM","JHN","JHO","JHP","JHQ","JHR","JHS","JHT","JHU","JHV",
"JHW","JHX","JHY","JHZ","JIA","JIB","JIC","JID","JIE","JIF","JIG","JIH","JII","JIJ","JIK","JIL",
"JIM","JIN","JIO","JIP","JIQ","JIR","JIS","JIT","JIU","JIV","JIW","JIX","JIY","JIZ","JJA","JJB",
"JJC","JJD","JJE","JJF","JJG","JJH","JJI","JJJ","JJK","JJL","JJM","JJN","JJO","JJP","JJQ","JJR",
"JJS","JJT","JJU","JJV","JJW","JJX","JJY","JJZ","JKA","JKB","JKC","JKD","JKE","JKF","JKG","JKH",
"JKI","JKJ","JKK","JKL","JKM","JKN","JKO","JKP","JKQ","JKR","JKS","JKT","JKU","JKV","JKW","JKX",
"JKY","JKZ","JLA","JLB","JLC","JLD","JLE","JLF","JLG","JLH","JLI","JLJ","JLK","JLL","JLM","JLN",
"JLO","JLP","JLQ","JLR","JLS","JLT","JLU","JLV","JLW","JLX","JLY","JLZ","JMA","JMB","JMC","JMD",
"JME","JMF","JMG","JMH","JMI","JMJ","JMK","JML","JMM","JMN","JMO","JMP","JMQ","JMR","JMS","JMT",
"JMU","JMV","JMW","JMX","JMY","JMZ","JNA","JNB","JNC","JND","JNE","JNF","JNG","JNH","JNI","JNJ",
"JNK","JNL","JNM","JNN","JNO","JNP","JNQ","JNR","JNS","JNT","JNU","JNV","JNW","JNX","JNY","JNZ",
"JOA","JOB","JOC","JOD","JOE","JOF","JOG","JOH","JOI","JOJ","JOK","JOL","JOM","JON","JOO","JOP",
"JOQ","JOR","JOS","JOT","JOU","JOV","JOW","JOX","JOY","JOZ","JPA","JPB","JPC","JPD","JPE","JPF",
"JPG","JPH","JPI","JPJ","JPK","JPL","JPM","JPN","JPO","JPP","JPQ","JPR","JPS","JPT","JPU","JPV",
"JPW","JPX","JPY","JPZ","JQA","JQB","JQC","JQD","JQE","JQF","JQG","JQH","JQI","JQJ","JQK","JQL",
"JQM","JQN","JQO","JQP","JQQ","JQR","JQS","JQT","JQU","JQV","JQW","JQX","JQY","JQZ","JRA","JRB",
"JRC","JRD","JRE","JRF","JRG","JRH","JRI","JRJ","JRK","JRL","JRM","JRN","JRO","JRP","JRQ","JRR",
"JRS","JRT","JRU","JRV","JRW","JRX","JRY","JRZ","JSA","JSB","JSC","JSD","JSE","JSF","JSG","JSH",
"JSI","JSJ","JSK","JSL","JSM","JSN","JSO","JSP","JSQ","JSR","JSS","JST","JSU","JSV","JSW","JSX",
"JSY","JSZ","JTA","JTB","JTC","JTD","JTE","JTF","JTG","JTH","JTI","JTJ","JTK","JTL","JTM","JTN",
"JTO","JTP","JTQ","JTR","JTS","JTT","JTU","JTV","JTW","JTX","JTY","JTZ","JUA","JUB","JUC","JUD",
"JUE","JUF","JUG","JUH","JUI","JUJ","JUK","JUL","JUM","JUN","JUO","JUP","JUQ","JUR","JUS","JUT",
"JUU","JUV","JUW","JUX","JUY","JUZ","JVA","JVB","JVC","JVD","JVE","JVF","JVG","JVH","JVI","JVJ",
"JVK","JVL","JVM","JVN","JVO","JVP","JVQ","JVR","JVS","JVT","JVU","JVV","JVW","JVX","JVY","JVZ",
"JWA","JWB","JWC","JWD","JWE","JWF","JWG","JWH","JWI","JWJ","JWK","JWL","JWM","JWN","JWO","JWP",
"JWQ","JWR","JWS","JWT","JWU","JWV","JWW","JWX","JWY","JWZ","JXA","JXB","JXC","JXD","JXE","JXF",
"JXG","JXH","JXI","JXJ","JXK","JXL","JXM","JXN","JXO","JXP","JXQ","JXR","JXS","JXT","JXU","JXV",
"JXW","JXX","JXY","JXZ","JYA","JYB","JYC","JYD","JYE","JYF","JYG","JYH","JYI","JYJ","JYK","JYL",
"JYM","JYN","JYO","JYP","JYQ","JYR","JYS","JYT","JYU","JYV","JYW","JYX","JYY","JYZ","JZA","JZB",
"JZC","JZD","JZE","JZF","JZG","JZH","JZI","JZJ","JZK","JZL","JZM","JZN","JZO","JZP","JZQ","JZR",
"JZS","JZT","JZU","JZV","JZW","JZX","JZY","JZZ","KAA","KAB","KAC","KAD","KAE","KAF","KAG","KAH",
"KAI","KAJ","KAK","KAL","KAM","KAN","KAO","KAP","KAQ","KAR","KAS","KAT","KAU","KAV","KAW","KAX",
"KAY","KAZ","KBA","KBB","KBC","KBD","KBE","KBF","KBG","KBH","KBI","KBJ","KBK","KBL","KBM","KBN",
"KBO","KBP","KBQ","KBR","KBS","KBT","KBU","KBV","KBW","KBX","KBY","KBZ","KCA","KCB","KCC","KCD",
"KCE","KCF","KCG","KCH","KCI","KCJ","KCK","KCL","KCM","KCN","KCO","KCP","KCQ","KCR","KCS","KCT",
"KCU","KCV","KCW","KCX","KCY","KCZ","KDA","KDB","KDC","KDD","KDE","KDF","KDG","KDH","KDI","KDJ",
"KDK","KDL","KDM","KDN","KDO","KDP","KDQ","KDR","KDS","KDT","KDU","KDV","KDW","KDX","KDY","KDZ",
"KEA","KEB","KEC","KED","KEE","KEF","KEG","KEH","KEI","KEJ","KEK","KEL","KEM","KEN","KEO","KEP",
"KEQ","KER","KES","KET","KEU","KEV","KEW","KEX","KEY","KEZ","KFA","KFB","KFC","KFD","KFE","KFF",
"KFG","KFH","KFI","KFJ","KFK","KFL","KFM","KFN","KFO","KFP","KFQ","KFR","KFS","KFT","KFU","KFV",
"KFW","KFX","KFY","KFZ","KGA","KGB","KGC","KGD","KGE","KGF","KGG","KGH","KGI","KGJ","KGK","KGL",
"KGM","KGN","KGO","KGP","KGQ","KGR","KGS","KGT","KGU","KGV","KGW","KGX","KGY","KGZ","KHA","KHB",
"KHC","KHD","KHE","KHF","KHG","KHH","KHI","KHJ","KHK","KHL","KHM","KHN","KHO","KHP","KHQ","KHR",
"KHS","KHT","KHU","KHV","KHW","KHX","KHY","KHZ","KIA","KIB","KIC","KID","KIE","KIF","KIG","KIH",
"KII","KIJ","KIK","KIL","KIM","KIN","KIO","KIP","KIQ","KIR","KIS","KIT","KIU","KIV","KIW","KIX",
"KIY","KIZ","KJA","KJB","KJC","KJD","KJE","KJF","KJG","KJH","KJI","KJJ","KJK","KJL","KJM","KJN",
"KJO","KJP","KJQ","KJR","KJS","KJT","KJU","KJV","KJW","KJX","KJY","KJZ","KKA","KKB","KKC","KKD",
"KKE","KKF","KKG","KKH","KKI","KKJ","KKK","KKL","KKM","KKN","KKO","KKP","KKQ","KKR","KKS","KKT",
"KKU","KKV","KKW","KKX","KKY","KKZ","KLA","KLB","KLC","KLD","KLE","KLF","KLG","KLH","KLI","KLJ",
"KLK","KLL","KLM","KLN","KLO","KLP","KLQ","KLR","KLS","KLT","KLU","KLV","KLW","KLX","KLY","KLZ",
"KMA","KMB","KMC","KMD","KME","KMF","KMG","KMH","KMI","KMJ","KMK","KML","KMM","KMN","KMO","KMP",
"KMQ","KMR","KMS","KMT","KMU","KMV","KMW","KMX","KMY","KMZ","KNA","KNB","KNC","KND","KNE","KNF",
"KNG","KNH","KNI","KNJ","KNK","KNL","KNM","KNN","KNO","KNP","KNQ","KNR","KNS","KNT","KNU","KNV",
"KNW","KNX","KNY","KNZ","KOA","KOB","KOC","KOD","KOE","KOF","KOG","KOH","KOI","KOJ","KOK","KOL",
"KOM","KON","KOO","KOP","KOQ","KOR","KOS","KOT","KOU","KOV","KOW","KOX","KOY","KOZ","KPA","KPB",
"KPC","KPD","KPE","KPF","KPG","KPH","KPI","KPJ","KPK","KPL","KPM","KPN","KPO","KPP","KPQ","KPR",
"KPS","KPT","KPU","KPV","KPW","KPX","KPY","KPZ","KQA","KQB","KQC","KQD","KQE","KQF","KQG","KQH",
"KQI","KQJ","KQK","KQL","KQM","KQN","KQO","KQP","KQQ","KQR","KQS","KQT","KQU","KQV","KQW","KQX",
"KQY","KQZ","KRA","KRB","KRC","KRD","KRE","KRF","KRG","KRH","KRI","KRJ","KRK","KRL","KRM","KRN",
"KRO","KRP","KRQ","KRR","KRS","KRT","KRU","KRV","KRW","KRX","KRY","KRZ","KSA","KSB","KSC","KSD",
"KSE","KSF","KSG","KSH","KSI","KSJ","KSK","KSL","KSM","KSN","KSO","KSP","KSQ","KSR","KSS","KST",
"KSU","KSV","KSW","KSX","KSY","KSZ","KTA","KTB","KTC","KTD","KTE","KTF","KTG","KTH","KTI","KTJ",
"KTK","KTL","KTM","KTN","KTO","KTP","KTQ","KTR","KTS","KTT","KTU","KTV","KTW","KTX","KTY","KTZ",
"KUA","KUB","KUC","KUD","KUE","KUF","KUG","KUH","KUI","KUJ","KUK","KUL","KUM","KUN","KUO","KUP",
"KUQ","KUR","KUS","KUT","KUU","KUV","KUW","KUX","KUY","KUZ","KVA","KVB","KVC","KVD","KVE","KVF",
"KVG","KVH","KVI","KVJ","KVK","KVL","KVM","KVN","KVO","KVP","KVQ","KVR","KVS","KVT","KVU","KVV",
"KVW","KVX","KVY","KVZ","KWA","KWB","KWC","KWD","KWE","KWF","KWG","KWH","KWI","KWJ","KWK","KWL",
"KWM","KWN","KWO","KWP","KWQ","KWR","KWS","KWT","KWU","KWV","KWW","KWX","KWY","KWZ","KXA","KXB",
"KXC","KXD","KXE","KXF","KXG","KXH","KXI","KXJ","KXK","KXL","KXM","KXN","KXO","KXP","KXQ","KXR",
"KXS","KXT","KXU","KXV","KXW","KXX","KXY","KXZ","KYA","KYB","KYC","KYD","KYE","KYF","KYG","KYH",
"KYI","KYJ","KYK","KYL","KYM","KYN","KYO","KYP","KYQ","KYR","KYS","KYT","KYU","KYV","KYW","KYX",
"KYY","KYZ","KZA","KZB","KZC","KZD","KZE","KZF","KZG","KZH","KZI","KZJ","KZK","KZL","KZM","KZN",
"KZO","KZP","KZQ","KZR","KZS","KZT","KZU","KZV","KZW","KZX","KZY","KZZ","LAA","LAB","LAC","LAD",
"LAE","LAF","LAG","LAH","LAI","LAJ","LAK","LAL","LAM","LAN","LAO","LAP","LAQ","LAR","LAS","LAT",
"LAU","LAV","LAW","LAX","LAY","LAZ","LBA","LBB","LBC","LBD","LBE","LBF","LBG","LBH","LBI","LBJ",
"LBK","LBL","LBM","LBN","LBO","LBP","LBQ","LBR","LBS","LBT","LBU","LBV","LBW","LBX","LBY","LBZ",
"LCA","LCB","LCC","LCD","LCE","LCF","LCG","LCH","LCI","LCJ","LCK","LCL","LCM","LCN","LCO","LCP",
"LCQ","LCR","LCS","LCT","LCU","LCV","LCW","LCX","LCY","LCZ","LDA","LDB","LDC","LDD","LDE","LDF",
"LDG","LDH","LDI","LDJ","LDK","LDL","LDM","LDN","LDO","LDP","LDQ","LDR","LDS","LDT","LDU","LDV",
"LDW","LDX","LDY","LDZ","LEA","LEB","LEC","LED","LEE","LEF","LEG","LEH","LEI","LEJ","LEK","LEL",
"LEM","LEN","LEO","LEP","LEQ","LER","LES","LET","LEU","LEV","LEW","LEX","LEY","LEZ","LFA","LFB",
"LFC","LFD","LFE","LFF","LFG","LFH","LFI","LFJ","LFK","LFL","LFM","LFN","LFO","LFP","LFQ","LFR",
"LFS","LFT","LFU","LFV","LFW","LFX","LFY","LFZ","LGA","LGB","LGC","LGD","LGE","LGF","LGG","LGH",
"LGI","LGJ","LGK","LGL","LGM","LGN","LGO","LGP","LGQ","LGR","LGS","LGT","LGU","LGV","LGW","LGX",
"LGY","LGZ","LHA","LHB","LHC","LHD","LHE","LHF","LHG","LHH","LHI","LHJ","LHK","LHL","LHM","LHN",
"LHO","LHP","LHQ","LHR","LHS","LHT","LHU","LHV","LHW","LHX","LHY","LHZ","LIA","LIB","LIC","LID",
"LIE","LIF","LIG","LIH","LII","LIJ","LIK","LIL","LIM","LIN","LIO","LIP","LIQ","LIR","LIS","LIT",
"LIU","LIV","LIW","LIX","LIY","LIZ","LJA","LJB","LJC","LJD","LJE","LJF","LJG","LJH","LJI","LJJ",
"LJK","LJL","LJM","LJN","LJO","LJP","LJQ","LJR","LJS","LJT","LJU","LJV","LJW","LJX","LJY","LJZ",
"LKA","LKB","LKC","LKD","LKE","LKF","LKG","LKH","LKI","LKJ","LKK","LKL","LKM","LKN","LKO","LKP",
"LKQ","LKR","LKS","LKT","LKU","LKV","LKW","LKX","LKY","LKZ","LLA","LLB","LLC","LLD","LLE","LLF",
"LLG","LLH","LLI","LLJ","LLK","LLL","LLM","LLN","LLO","LLP","LLQ","LLR","LLS","LLT","LLU","LLV",
"LLW","LLX","LLY","LLZ","LMA","LMB","LMC","LMD","LME","LMF","LMG","LMH","LMI","LMJ","LMK","LML",
"LMM","LMN","LMO","LMP","LMQ","LMR","LMS","LMT","LMU","LMV","LMW","LMX","LMY","LMZ","LNA","LNB",
"LNC","LND","LNE","LNF","LNG","LNH","LNI","LNJ","LNK","LNL","LNM","LNN","LNO","LNP","LNQ","LNR",
"LNS","LNT","LNU","LNV","LNW","LNX","LNY","LNZ","LOA","LOB","LOC","LOD","LOE","LOF","LOG","LOH",
"LOI","LOJ","LOK","LOL","LOM","LON","LOO","LOP","LOQ","LOR","LOS","LOT","LOU","LOV","LOW","LOX",
"LOY","LOZ","LPA","LPB","LPC","LPD","LPE","LPF","LPG","LPH","LPI","LPJ","LPK","LPL","LPM","LPN",
"LPO","LPP","LPQ","LPR","LPS","LPT","LPU","LPV","LPW","LPX","LPY","LPZ","LQA","LQB","LQC","LQD",
"LQE","LQF","LQG","LQH","LQI","LQJ","LQK","LQL","LQM","LQN","LQO","LQP","LQQ","LQR","LQS","LQT",
"LQU","LQV","LQW","LQX","LQY","LQZ","LRA","LRB","LRC","LRD","LRE","LRF","LRG","LRH","LRI","LRJ",
"LRK","LRL","LRM","LRN","LRO","LRP","LRQ","LRR","LRS","LRT","LRU","LRV","LRW","LRX","LRY","LRZ",
"LSA","LSB","LSC","LSD","LSE","LSF","LSG","LSH","LSI","LSJ","LSK","LSL","LSM","LSN","LSO","LSP",
"LSQ","LSR","LSS","LST","LSU","LSV","LSW","LSX","LSY","LSZ","LTA","LTB","LTC","LTD","LTE","LTF",
"LTG","LTH","LTI","LTJ","LTK","LTL","LTM","LTN","LTO","LTP","LTQ","LTR","LTS","LTT","LTU","LTV",
"LTW","LTX","LTY","LTZ","LUA","LUB","LUC","LUD","LUE","LUF","LUG","LUH","LUI","LUJ","LUK","LUL",
"LUM","LUN","LUO","LUP","LUQ","LUR","LUS","LUT","LUU","LUV","LUW","LUX","LUY","LUZ","LVA","LVB",
"LVC","LVD","LVE","LVF","LVG","LVH","LVI","LVJ","LVK","LVL","LVM","LVN","LVO","LVP","LVQ","LVR",
"LVS","LVT","LVU","LVV","LVW","LVX","LVY","LVZ","LWA","LWB","LWC","LWD","LWE","LWF","LWG","LWH",
"LWI","LWJ","LWK","LWL","LWM","LWN","LWO","LWP","LWQ","LWR","LWS","LWT","LWU","LWV","LWW","LWX",
"LWY","LWZ","LXA","LXB","LXC","LXD","LXE","LXF","LXG","LXH","LXI","LXJ","LXK","LXL","LXM","LXN",
"LXO","LXP","LXQ","LXR","LXS","LXT","LXU","LXV","LXW","LXX","LXY","LXZ","LYA","LYB","LYC","LYD",
"LYE","LYF","LYG","LYH","LYI","LYJ","LYK","LYL","LYM","LYN","LYO","LYP","LYQ","LYR","LYS","LYT",
"LYU","LYV","LYW","LYX","LYY","LYZ","LZA","LZB","LZC","LZD","LZE","LZF","LZG","LZH","LZI","LZJ",
"LZK","LZL","LZM","LZN","LZO","LZP","LZQ","LZR","LZS","LZT","LZU","LZV","LZW","LZX","LZY","LZZ",
"MAA","MAB","MAC","MAD","MAE","MAF","MAG","MAH","MAI","MAJ","MAK","MAL","MAM","MAN","MAO","MAP",
"MAQ","MAR","MAS","MAT","MAU","MAV","MAW","MAX","MAY","MAZ","MBA","MBB","MBC","MBD","MBE","MBF",
"MBG","MBH","MBI","MBJ","MBK","MBL","MBM","MBN","MBO","MBP","MBQ","MBR","MBS","MBT","MBU","MBV",
"MBW","MBX","MBY","MBZ","MCA","MCB","MCC","MCD","MCE","MCF","MCG","MCH","MCI","MCJ","MCK","MCL",
"MCM","MCN","MCO","MCP","MCQ","MCR","MCS","MCT","MCU","MCV","MCW","MCX","MCY","MCZ","MDA","MDB",
"MDC","MDD","MDE","MDF","MDG","MDH","MDI","MDJ","MDK","MDL","MDM","MDN","MDO","MDP","MDQ","MDR",
"MDS","MDT","MDU","MDV","MDW","MDX","MDY","MDZ","MEA","MEB","MEC","MED","MEE","MEF","MEG","MEH",
"MEI","MEJ","MEK","MEL","MEM","MEN","MEO","MEP","MEQ","MER","MES","MET","MEU","MEV","MEW","MEX",
"MEY","MEZ","MFA","MFB","MFC","MFD","MFE","MFF","MFG","MFH","MFI","MFJ","MFK","MFL","MFM","MFN",
"MFO","MFP","MFQ","MFR","MFS","MFT","MFU","MFV","MFW","MFX","MFY","MFZ","MGA","MGB","MGC","MGD",
"MGE","MGF","MGG","MGH","MGI","MGJ","MGK","MGL","MGM","MGN","MGO","MGP","MGQ","MGR","MGS","MGT",
"MGU","MGV","MGW","MGX","MGY","MGZ","MHA","MHB","MHC","MHD","MHE","MHF","MHG","MHH","MHI","MHJ",
"MHK","MHL","MHM","MHN","MHO","MHP","MHQ","MHR","MHS","MHT","MHU","MHV","MHW","MHX","MHY","MHZ",
"MIA","MIB","MIC","MID","MIE","MIF","MIG","MIH","MII","MIJ","MIK","MIL","MIM","MIN","MIO","MIP",
"MIQ","MIR","MIS","MIT","MIU","MIV","MIW","MIX","MIY","MIZ","MJA","MJB","MJC","MJD","MJE","MJF",
"MJG","MJH","MJI","MJJ","MJK","MJL","MJM","MJN","MJO","MJP","MJQ","MJR","MJS","MJT","MJU","MJV",
"MJW","MJX","MJY","MJZ","MKA","MKB","MKC","MKD","MKE","MKF","MKG","MKH","MKI","MKJ","MKK","MKL",
"MKM","MKN","MKO","MKP","MKQ","MKR","MKS","MKT","MKU","MKV","MKW","MKX","MKY","MKZ","MLA","MLB",
"MLC","MLD","MLE","MLF","MLG","MLH","MLI","MLJ","MLK","MLL","MLM","MLN","MLO","MLP","MLQ","MLR",
"MLS","MLT","MLU","MLV","MLW","MLX","MLY","MLZ","MMA","MMB","MMC","MMD","MME","MMF","MMG","MMH",
"MMI","MMJ","MMK","MML","MMM","MMN","MMO","MMP","MMQ","MMR","MMS","MMT","MMU","MMV","MMW","MMX",
"MMY","MMZ","MNA","MNB","MNC","MND","MNE","MNF","MNG","MNH","MNI","MNJ","MNK","MNL","MNM","MNN",
"MNO","MNP","MNQ","MNR","MNS","MNT","MNU","MNV","MNW","MNX","MNY","MNZ","MOA","MOB","MOC","MOD",
"MOE","MOF","MOG","MOH","MOI","MOJ","MOK","MOL","MOM","MON","MOO","MOP","MOQ","MOR","MOS","MOT",
"MOU","MOV","MOW","MOX","MOY","MOZ","MPA","MPB","MPC","MPD","MPE","MPF","MPG","MPH","MPI","MPJ",
"MPK","MPL","MPM","MPN","MPO","MPP","MPQ","MPR","MPS","MPT","MPU","MPV","MPW","MPX","MPY","MPZ",
"MQA","MQB","MQC","MQD","MQE","MQF","MQG","MQH","MQI","MQJ","MQK","MQL","MQM","MQN","MQO","MQP",
"MQQ","MQR","MQS","MQT","MQU","MQV","MQW","MQX","MQY","MQZ","MRA","MRB","MRC","MRD","MRE","MRF",
"MRG","MRH","MRI","MRJ","MRK","MRL","MRM","MRN","MRO","MRP","MRQ","MRR","MRS","MRT","MRU","MRV",
"MRW","MRX","MRY","MRZ","MSA","MSB","MSC","MSD","MSE","MSF","MSG","MSH","MSI","MSJ","MSK","MSL",
"MSM","MSN","MSO","MSP","MSQ","MSR","MSS","MST","MSU","MSV","MSW","MSX","MSY","MSZ","MTA","MTB",
"MTC","MTD","MTE","MTF","MTG","MTH","MTI","MTJ","MTK","MTL","MTM","MTN","MTO","MTP","MTQ","MTR",
"MTS","MTT","MTU","MTV","MTW","MTX","MTY","MTZ","MUA","MUB","MUC","MUD","MUE","MUF","MUG","MUH",
"MUI","MUJ","MUK","MUL","MUM","MUN","MUO","MUP","MUQ","MUR","MUS","MUT","MUU","MUV","MUW","MUX",
"MUY","MUZ","MVA","MVB","MVC","MVD","MVE","MVF","MVG","MVH","MVI","MVJ","MVK","MVL","MVM","MVN",
"MVO","MVP","MVQ","MVR","MVS","MVT","MVU","MVV","MVW","MVX","MVY","MVZ","MWA","MWB","MWC","MWD",
"MWE","MWF","MWG","MWH","MWI","MWJ","MWK","MWL","MWM","MWN","MWO","MWP","MWQ","MWR","MWS","MWT",
"MWU","MWV","MWW","MWX","MWY","MWZ","MXA","MXB","MXC","MXD","MXE","MXF","MXG","MXH","MXI","MXJ",
"MXK","MXL","MXM","MXN","MXO","MXP","MXQ","MXR","MXS","MXT","MXU","MXV","MXW","MXX","MXY","MXZ",
"MYA","MYB","MYC","MYD","MYE","MYF","MYG","MYH","MYI","MYJ","MYK","MYL","MYM","MYN","MYO","MYP",
"MYQ","MYR","MYS","MYT","MYU","MYV","MYW","MYX","MYY","MYZ","MZA","MZB","MZC","MZD","MZE","MZF",
"MZG","MZH","MZI","MZJ","MZK","MZL","MZM","MZN","MZO","MZP","MZQ","MZR","MZS","MZT","MZU","MZV",
"MZW","MZX","MZY","MZZ","NAA","NAB","NAC","NAD","NAE","NAF","NAG","NAH","NAI","NAJ","NAK","NAL",
"NAM","NAN","NAO","NAP","NAQ","NAR","NAS","NAT","NAU","NAV","NAW","NAX","NAY","NAZ","NBA","NBB",
"NBC","NBD","NBE","NBF","NBG","NBH","NBI","NBJ","NBK","NBL","NBM","NBN","NBO","NBP","NBQ","NBR",
"NBS","NBT","NBU","NBV","NBW","NBX","NBY","NBZ","NCA","NCB","NCC","NCD","NCE","NCF","NCG","NCH",
"NCI","NCJ","NCK","NCL","NCM","NCN","NCO","NCP","NCQ","NCR","NCS","NCT","NCU","NCV","NCW","NCX",
"NCY","NCZ","NDA","NDB","NDC","NDD","NDE","NDF","NDG","NDH","NDI","NDJ","NDK","NDL","NDM","NDN",
"NDO","NDP","NDQ","NDR","NDS","NDT","NDU","NDV","NDW","NDX","NDY","NDZ","NEA","NEB","NEC","NED",
"NEE","NEF","NEG","NEH","NEI","NEJ","NEK","NEL","NEM","NEN","NEO","NEP","NEQ","NER","NES","NET",
"NEU","NEV","NEW","NEX","NEY","NEZ","NFA","NFB","NFC","NFD","NFE","NFF","NFG","NFH","NFI","NFJ",
"NFK","NFL","NFM","NFN","NFO","NFP","NFQ","NFR","NFS","NFT","NFU","NFV","NFW","NFX","NFY","NFZ",
"NGA","NGB","NGC","NGD","NGE","NGF","NGG","NGH","NGI","NGJ","NGK","NGL","NGM","NGN","NGO","NGP",
"NGQ","NGR","NGS","NGT","NGU","NGV","NGW","NGX","NGY","NGZ","NHA","NHB","NHC","NHD","NHE","NHF",
"NHG","NHH","NHI","NHJ","NHK","NHL","NHM","NHN","NHO","NHP","NHQ","NHR","NHS","NHT","NHU","NHV",
"NHW","NHX","NHY","NHZ","NIA","NIB","NIC","NID","NIE","NIF","NIG","NIH","NII","NIJ","NIK","NIL",
"NIM","NIN","NIO","NIP","NIQ","NIR","NIS","NIT","NIU","NIV","NIW","NIX","NIY","NIZ","NJA","NJB",
"NJC","NJD","NJE","NJF","NJG","NJH","NJI","NJJ","NJK","NJL","NJM","NJN","NJO","NJP","NJQ","NJR",
"NJS","NJT","NJU","NJV","NJW","NJX","NJY","NJZ","NKA","NKB","NKC","NKD","NKE","NKF","NKG","NKH",
"NKI","NKJ","NKK","NKL","NKM","NKN","NKO","NKP","NKQ","NKR","NKS","NKT","NKU","NKV","NKW","NKX",
"NKY","NKZ","NLA","NLB","NLC","NLD","NLE","NLF","NLG","NLH","NLI","NLJ","NLK","NLL","NLM","NLN",
"NLO","NLP","NLQ","NLR","NLS","NLT","NLU","NLV","NLW","NLX","NLY","NLZ","NMA","NMB","NMC","NMD",
"NME","NMF","NMG","NMH","NMI","NMJ","NMK","NML","NMM","NMN","NMO","NMP","NMQ","NMR","NMS","NMT",
"NMU","NMV","NMW","NMX","NMY","NMZ","NNA","NNB","NNC","NND","NNE","NNF","NNG","NNH","NNI","NNJ",
"NNK","NNL","NNM","NNN","NNO","NNP","NNQ","NNR","NNS","NNT","NNU","NNV","NNW","NNX","NNY","NNZ",
"NOA","NOB","NOC","NOD","NOE","NOF","NOG","NOH","NOI","NOJ","NOK","NOL","NOM","NON","NOO","NOP",
"NOQ","NOR","NOS","NOT","NOU","NOV","NOW","NOX","NOY","NOZ","NPA","NPB","NPC","NPD","NPE","NPF",
"NPG","NPH","NPI","NPJ","NPK","NPL","NPM","NPN","NPO","NPP","NPQ","NPR","NPS","NPT","NPU","NPV",
"NPW","NPX","NPY","NPZ","NQA","NQB","NQC","NQD","NQE","NQF","NQG","NQH","NQI","NQJ","NQK","NQL",
"NQM","NQN","NQO","NQP","NQQ","NQR","NQS","NQT","NQU","NQV","NQW","NQX","NQY","NQZ","NRA","NRB",
"NRC","NRD","NRE","NRF","NRG","NRH","NRI","NRJ","NRK","NRL","NRM","NRN","NRO","NRP","NRQ","NRR",
"NRS","NRT","NRU","NRV","NRW","NRX","NRY","NRZ","NSA","NSB","NSC","NSD","NSE","NSF","NSG","NSH",
"NSI","NSJ","NSK","NSL","NSM","NSN","NSO","NSP","NSQ","NSR","NSS","NST","NSU","NSV","NSW","NSX",
"NSY","NSZ","NTA","NTB","NTC","NTD","NTE","NTF","NTG","NTH","NTI","NTJ","NTK","NTL","NTM","NTN",
"NTO","NTP","NTQ","NTR","NTS","NTT","NTU","NTV","NTW","NTX","NTY","NTZ","NUA","NUB","NUC","NUD",
"NUE","NUF","NUG","NUH","NUI","NUJ","NUK","NUL","NUM","NUN","NUO","NUP","NUQ","NUR","NUS","NUT",
"NUU","NUV","NUW","NUX","NUY","NUZ","NVA","NVB","NVC","NVD","NVE","NVF","NVG","NVH","NVI","NVJ",
"NVK","NVL","NVM","NVN","NVO","NVP","NVQ","NVR","NVS","NVT","NVU","NVV","NVW","NVX","NVY","NVZ",
"NWA","NWB","NWC","NWD","NWE","NWF","NWG","NWH","NWI","NWJ","NWK","NWL","NWM","NWN","NWO","NWP",
"NWQ","NWR","NWS","NWT","NWU","NWV","NWW","NWX","NWY","NWZ","NXA","NXB","NXC","NXD","NXE","NXF",
"NXG","NXH","NXI","NXJ","NXK","NXL","NXM","NXN","NXO","NXP","NXQ","NXR","NXS","NXT","NXU","NXV",
"NXW","NXX","NXY","NXZ","NYA","NYB","NYC","NYD","NYE","NYF","NYG","NYH","NYI","NYJ","NYK","NYL",
"NYM","NYN","NYO","NYP","NYQ","NYR","NYS","NYT","NYU","NYV","NYW","NYX","NYY","NYZ","NZA","NZB",
"NZC","NZD","NZE","NZF","NZG","NZH","NZI","NZJ","NZK","NZL","NZM","NZN","NZO","NZP","NZQ","NZR",
"NZS","NZT","NZU","NZV","NZW","NZX","NZY","NZZ","OAA","OAB","OAC","OAD","OAE","OAF","OAG","OAH",
"OAI","OAJ","OAK","OAL","OAM","OAN","OAO","OAP","OAQ","OAR","OAS","OAT","OAU","OAV","OAW","OAX",
"OAY","OAZ","OBA","OBB","OBC","OBD","OBE","OBF","OBG","OBH","OBI","OBJ","OBK","OBL","OBM","OBN",
"OBO","OBP","OBQ","OBR","OBS","OBT","OBU","OBV","OBW","OBX","OBY","OBZ","OCA","OCB","OCC","OCD",
"OCE","OCF","OCG","OCH","OCI","OCJ","OCK","OCL","OCM","OCN","OCO","OCP","OCQ","OCR","OCS","OCT",
"OCU","OCV","OCW","OCX","OCY","OCZ","ODA","ODB","ODC","ODD","ODE","ODF","ODG","ODH","ODI","ODJ",
"ODK","ODL","ODM","ODN","ODO","ODP","ODQ","ODR","ODS","ODT","ODU","ODV","ODW","ODX","ODY","ODZ",
"OEA","OEB","OEC","OED","OEE","OEF","OEG","OEH","OEI","OEJ","OEK","OEL","OEM","OEN","OEO","OEP",
"OEQ","OER","OES","OET","OEU","OEV","OEW","OEX","OEY","OEZ","OFA","OFB","OFC","OFD","OFE","OFF",
"OFG","OFH","OFI","OFJ","OFK","OFL","OFM","OFN","OFO","OFP","OFQ","OFR","OFS","OFT","OFU","OFV",
"OFW","OFX","OFY","OFZ","OGA","OGB","OGC","OGD","OGE","OGF","OGG","OGH","OGI","OGJ","OGK","OGL",
"OGM","OGN","OGO","OGP","OGQ","OGR","OGS","OGT","OGU","OGV","OGW","OGX","OGY","OGZ","OHA","OHB",
"OHC","OHD","OHE","OHF","OHG","OHH","OHI","OHJ","OHK","OHL","OHM","OHN","OHO","OHP","OHQ","OHR",
"OHS","OHT","OHU","OHV","OHW","OHX","OHY","OHZ","OIA","OIB","OIC","OID","OIE","OIF","OIG","OIH",
"OII","OIJ","OIK","OIL","OIM","OIN","OIO","OIP","OIQ","OIR","OIS","OIT","OIU","OIV","OIW","OIX",
"OIY","OIZ","OJA","OJB","OJC","OJD","OJE","OJF","OJG","OJH","OJI","OJJ","OJK","OJL","OJM","OJN",
"OJO","OJP","OJQ","OJR","OJS","OJT","OJU","OJV","OJW","OJX","OJY","OJZ","OKA","OKB","OKC","OKD",
"OKE","OKF","OKG","OKH","OKI","OKJ","OKK","OKL","OKM","OKN","OKO","OKP","OKQ","OKR","OKS","OKT",
"OKU","OKV","OKW","OKX","OKY","OKZ","OLA","OLB","OLC","OLD","OLE","OLF","OLG","OLH","OLI","OLJ",
"OLK","OLL","OLM","OLN","OLO","OLP","OLQ","OLR","OLS","OLT","OLU","OLV","OLW","OLX","OLY","OLZ",
"OMA","OMB","OMC","OMD","OME","OMF","OMG","OMH","OMI","OMJ","OMK","OML","OMM","OMN","OMO","OMP",
"OMQ","OMR","OMS","OMT","OMU","OMV","OMW","OMX","OMY","OMZ","ONA","ONB","ONC","OND","ONE","ONF",
"ONG","ONH","ONI","ONJ","ONK","ONL","ONM","ONN","ONO","ONP","ONQ","ONR","ONS","ONT","ONU","ONV",
"ONW","ONX","ONY","ONZ","OOA","OOB","OOC","OOD","OOE","OOF","OOG","OOH","OOI","OOJ","OOK","OOL",
"OOM","OON","OOO","OOP","OOQ","OOR","OOS","OOT","OOU","OOV","OOW","OOX","OOY","OOZ","OPA","OPB",
"OPC","OPD","OPE","OPF","OPG","OPH","OPI","OPJ","OPK","OPL","OPM","OPN","OPO","OPP","OPQ","OPR",
"OPS","OPT","OPU","OPV","OPW","OPX","OPY","OPZ","OQA","OQB","OQC","OQD","OQE","OQF","OQG","OQH",
"OQI","OQJ","OQK","OQL","OQM","OQN","OQO","OQP","OQQ","OQR","OQS","OQT","OQU","OQV","OQW","OQX",
"OQY","OQZ","ORA","ORB","ORC","ORD","ORE","ORF","ORG","ORH","ORI","ORJ","ORK","ORL","ORM","ORN",
"ORO","ORP","ORQ","ORR","ORS","ORT","ORU","ORV","ORW","ORX","ORY","ORZ","OSA","OSB","OSC","OSD",
"OSE","OSF","OSG","OSH","OSI","OSJ","OSK","OSL","OSM","OSN","OSO","OSP","OSQ","OSR","OSS","OST",
"OSU","OSV","OSW","OSX","OSY","OSZ","OTA","OTB","OTC","OTD","OTE","OTF","OTG","OTH","OTI","OTJ",
"OTK","OTL","OTM","OTN","OTO","OTP","OTQ","OTR","OTS","OTT","OTU","OTV","OTW","OTX","OTY","OTZ",
"OUA","OUB","OUC","OUD","OUE","OUF","OUG","OUH","OUI","OUJ","OUK","OUL","OUM","OUN","OUO","OUP",
"OUQ","OUR","OUS","OUT","OUU","OUV","OUW","OUX","OUY","OUZ","OVA","OVB","OVC","OVD","OVE","OVF",
"OVG","OVH","OVI","OVJ","OVK","OVL","OVM","OVN","OVO","OVP","OVQ","OVR","OVS","OVT","OVU","OVV",
"OVW","OVX","OVY","OVZ","OWA","OWB","OWC","OWD","OWE","OWF","OWG","OWH","OWI","OWJ","OWK","OWL",
"OWM","OWN","OWO","OWP","OWQ","OWR","OWS","OWT","OWU","OWV","OWW","OWX","OWY","OWZ","OXA","OXB",
"OXC","OXD","OXE","OXF","OXG","OXH","OXI","OXJ","OXK","OXL","OXM","OXN","OXO","OXP","OXQ","OXR",
"OXS","OXT","OXU","OXV","OXW","OXX","OXY","OXZ","OYA","OYB","OYC","OYD","OYE","OYF","OYG","OYH",
"OYI","OYJ","OYK","OYL","OYM","OYN","OYO","OYP","OYQ","OYR","OYS","OYT","OYU","OYV","OYW","OYX",
"OYY","OYZ","OZA","OZB","OZC","OZD","OZE","OZF","OZG","OZH","OZI","OZJ","OZK","OZL","OZM","OZN",
"OZO","OZP","OZQ","OZR","OZS","OZT","OZU","OZV","OZW","OZX","OZY","OZZ","PAA","PAB","PAC","PAD",
"PAE","PAF","PAG","PAH","PAI","PAJ","PAK","PAL","PAM","PAN","PAO","PAP","PAQ","PAR","PAS","PAT",
"PAU","PAV","PAW","PAX","PAY","PAZ","PBA","PBB","PBC","PBD","PBE","PBF","PBG","PBH","PBI","PBJ",
"PBK","PBL","PBM","PBN","PBO","PBP","PBQ","PBR","PBS","PBT","PBU","PBV","PBW","PBX","PBY","PBZ",
"PCA","PCB","PCC","PCD","PCE","PCF","PCG","PCH","PCI","PCJ","PCK","PCL","PCM","PCN","PCO","PCP",
"PCQ","PCR","PCS","PCT","PCU","PCV","PCW","PCX","PCY","PCZ","PDA","PDB","PDC","PDD","PDE","PDF",
"PDG","PDH","PDI","PDJ","PDK","PDL","PDM","PDN","PDO","PDP","PDQ","PDR","PDS","PDT","PDU","PDV",
"PDW","PDX","PDY","PDZ","PEA","PEB","PEC","PED","PEE","PEF","PEG","PEH","PEI","PEJ","PEK","PEL",
"PEM","PEN","PEO","PEP","PEQ","PER","PES","PET","PEU","PEV","PEW","PEX","PEY","PEZ","PFA","PFB",
"PFC","PFD","PFE","PFF","PFG","PFH","PFI","PFJ","PFK","PFL","PFM","PFN","PFO","PFP","PFQ","PFR",
"PFS","PFT","PFU","PFV","PFW","PFX","PFY","PFZ","PGA","PGB","PGC","PGD","PGE","PGF","PGG","PGH",
"PGI","PGJ","PGK","PGL","PGM","PGN","PGO","PGP","PGQ","PGR","PGS","PGT","PGU","PGV","PGW","PGX",
"PGY","PGZ","PHA","PHB","PHC","PHD","PHE","PHF","PHG","PHH","PHI","PHJ","PHK","PHL","PHM","PHN",
"PHO","PHP","PHQ","PHR","PHS","PHT","PHU","PHV","PHW","PHX","PHY","PHZ","PIA","PIB","PIC","PID",
"PIE","PIF","PIG","PIH","PII","PIJ","PIK","PIL","PIM","PIN","PIO","PIP","PIQ","PIR","PIS","PIT",
"PIU","PIV","PIW","PIX","PIY","PIZ","PJA","PJB","PJC","PJD","PJE","PJF","PJG","PJH","PJI","PJJ",
"PJK","PJL","PJM","PJN","PJO","PJP","PJQ","PJR","PJS","PJT","PJU","PJV","PJW","PJX","PJY","PJZ",
"PKA","PKB","PKC","PKD","PKE","PKF","PKG","PKH","PKI","PKJ","PKK","PKL","PKM","PKN","PKO","PKP",
"PKQ","PKR","PKS","PKT","PKU","PKV","PKW","PKX","PKY","PKZ","PLA","PLB","PLC","PLD","PLE","PLF",
"PLG","PLH","PLI","PLJ","PLK","PLL","PLM","PLN","PLO","PLP","PLQ","PLR","PLS","PLT","PLU","PLV",
"PLW","PLX","PLY","PLZ","PMA","PMB","PMC","PMD","PME","PMF","PMG","PMH","PMI","PMJ","PMK","PML",
"PMM","PMN","PMO","PMP","PMQ","PMR","PMS","PMT","PMU","PMV","PMW","PMX","PMY","PMZ","PNA","PNB",
"PNC","PND","PNE","PNF","PNG","PNH","PNI","PNJ","PNK","PNL","PNM","PNN","PNO","PNP","PNQ","PNR",
"PNS","PNT","PNU","PNV","PNW","PNX","PNY","PNZ","POA","POB","POC","POD","POE","POF","POG","POH",
"POI","POJ","POK","POL","POM","PON","POO","POP","POQ","POR","POS","POT","POU","POV","POW","POX",
"POY","POZ","PPA","PPB","PPC","PPD","PPE","PPF","PPG","PPH","PPI","PPJ","PPK","PPL","PPM","PPN",
"PPO","PPP","PPQ","PPR","PPS","PPT","PPU","PPV","PPW","PPX","PPY","PPZ","PQA","PQB","PQC","PQD",
"PQE","PQF","PQG","PQH","PQI","PQJ","PQK","PQL","PQM","PQN","PQO","PQP","PQQ","PQR","PQS","PQT",
"PQU","PQV","PQW","PQX","PQY","PQZ","PRA","PRB","PRC","PRD","PRE","PRF","PRG","PRH","PRI","PRJ",
"PRK","PRL","PRM","PRN","PRO","PRP","PRQ","PRR","PRS","PRT","PRU","PRV","PRW","PRX","PRY","PRZ",
"PSA","PSB","PSC","PSD","PSE","PSF","PSG","PSH","PSI","PSJ","PSK","PSL","PSM","PSN","PSO","PSP",
"PSQ","PSR","PSS","PST","PSU","PSV","PSW","PSX","PSY","PSZ","PTA","PTB","PTC","PTD","PTE","PTF",
"PTG","PTH","PTI","PTJ","PTK","PTL","PTM","PTN","PTO","PTP","PTQ","PTR","PTS","PTT","PTU","PTV",
"PTW","PTX","PTY","PTZ","PUA","PUB","PUC","PUD","PUE","PUF","PUG","PUH","PUI","PUJ","PUK","PUL",
"PUM","PUN","PUO","PUP","PUQ","PUR","PUS","PUT","PUU","PUV","PUW","PUX","PUY","PUZ","PVA","PVB",
"PVC","PVD","PVE","PVF","PVG","PVH","PVI","PVJ","PVK","PVL","PVM","PVN","PVO","PVP","PVQ","PVR",
"PVS","PVT","PVU","PVV","PVW","PVX","PVY","PVZ","PWA","PWB","PWC","PWD","PWE","PWF","PWG","PWH",
"PWI","PWJ","PWK","PWL","PWM","PWN","PWO","PWP","PWQ","PWR","PWS","PWT","PWU","PWV","PWW","PWX",
"PWY","PWZ","PXA","PXB","PXC","PXD","PXE","PXF","PXG","PXH","PXI","PXJ","PXK","PXL","PXM","PXN",
"PXO","PXP","PXQ","PXR","PXS","PXT","PXU","PXV","PXW","PXX","PXY","PXZ","PYA","PYB","PYC","PYD",
"PYE","PYF","PYG","PYH","PYI","PYJ","PYK","PYL","PYM","PYN","PYO","PYP","PYQ","PYR","PYS","PYT",
"PYU","PYV","PYW","PYX","PYY","PYZ","PZA","PZB","PZC","PZD","PZE","PZF","PZG","PZH","PZI","PZJ",
"PZK","PZL","PZM","PZN","PZO","PZP","PZQ","PZR","PZS","PZT","PZU","PZV","PZW","PZX","PZY","PZZ",
"QAA","QAB","QAC","QAD","QAE","QAF","QAG","QAH","QAI","QAJ","QAK","QAL","QAM","QAN","QAO","QAP",
"QAQ","QAR","QAS","QAT","QAU","QAV","QAW","QAX","QAY","QAZ","QBA","QBB","QBC","QBD","QBE","QBF",
"QBG","QBH","QBI","QBJ","QBK","QBL","QBM","QBN","QBO","QBP","QBQ","QBR","QBS","QBT","QBU","QBV",
"QBW","QBX","QBY","QBZ","QCA","QCB","QCC","QCD","QCE","QCF","QCG","QCH","QCI","QCJ","QCK","QCL",
"QCM","QCN","QCO","QCP","QCQ","QCR","QCS","QCT","QCU","QCV","QCW","QCX","QCY","QCZ","QDA","QDB",
"QDC","QDD","QDE","QDF","QDG","QDH","QDI","QDJ","QDK","QDL","QDM","QDN","QDO","QDP","QDQ","QDR",
"QDS","QDT","QDU","QDV","QDW","QDX","QDY","QDZ","QEA","QEB","QEC","QED","QEE","QEF","QEG","QEH",
"QEI","QEJ","QEK","QEL","QEM","QEN","QEO","QEP","QEQ","QER","QES","QET","QEU","QEV","QEW","QEX",
"QEY","QEZ","QFA","QFB","QFC","QFD","QFE","QFF","QFG","QFH","QFI","QFJ","QFK","QFL","QFM","QFN",
"QFO","QFP","QFQ","QFR","QFS","QFT","QFU","QFV","QFW","QFX","QFY","QFZ","QGA","QGB","QGC","QGD",
"QGE","QGF","QGG","QGH","QGI","QGJ","QGK","QGL","QGM","QGN","QGO","QGP","QGQ","QGR","QGS","QGT",
"QGU","QGV","QGW","QGX","QGY","QGZ","QHA","QHB","QHC","QHD","QHE","QHF","QHG","QHH","QHI","QHJ",
"QHK","QHL","QHM","QHN","QHO","QHP","QHQ","QHR","QHS","QHT","QHU","QHV","QHW","QHX","QHY","QHZ",
"QIA","QIB","QIC","QID","QIE","QIF","QIG","QIH","QII","QIJ","QIK","QIL","QIM","QIN","QIO","QIP",
"QIQ","QIR","QIS","QIT","QIU","QIV","QIW","QIX","QIY","QIZ","QJA","QJB","QJC","QJD","QJE","QJF",
"QJG","QJH","QJI","QJJ","QJK","QJL","QJM","QJN","QJO","QJP","QJQ","QJR","QJS","QJT","QJU","QJV",
"QJW","QJX","QJY","QJZ","QKA","QKB","QKC","QKD","QKE","QKF","QKG","QKH","QKI","QKJ","QKK","QKL",
"QKM","QKN","QKO","QKP","QKQ","QKR","QKS","QKT","QKU","QKV","QKW","QKX","QKY","QKZ","QLA","QLB",
"QLC","QLD","QLE","QLF","QLG","QLH","QLI","QLJ","QLK","QLL","QLM","QLN","QLO","QLP","QLQ","QLR",
"QLS","QLT","QLU","QLV","QLW","QLX","QLY","QLZ","QMA","QMB","QMC","QMD","QME","QMF","QMG","QMH",
"QMI","QMJ","QMK","QML","QMM","QMN","QMO","QMP","QMQ","QMR","QMS","QMT","QMU","QMV","QMW","QMX",
"QMY","QMZ","QNA","QNB","QNC","QND","QNE","QNF","QNG","QNH","QNI","QNJ","QNK","QNL","QNM","QNN",
"QNO","QNP","QNQ","QNR","QNS","QNT","QNU","QNV","QNW","QNX","QNY","QNZ","QOA","QOB","QOC","QOD",
"QOE","QOF","QOG","QOH","QOI","QOJ","QOK","QOL","QOM","QON","QOO","QOP","QOQ","QOR","QOS","QOT",
"QOU","QOV","QOW","QOX","QOY","QOZ","QPA","QPB","QPC","QPD","QPE","QPF","QPG","QPH","QPI","QPJ",
"QPK","QPL","QPM","QPN","QPO","QPP","QPQ","QPR","QPS","QPT","QPU","QPV","QPW","QPX","QPY","QPZ",
"QQA","QQB","QQC","QQD","QQE","QQF","QQG","QQH","QQI","QQJ","QQK","QQL","QQM","QQN","QQO","QQP",
"QQQ","QQR","QQS","QQT","QQU","QQV","QQW","QQX","QQY","QQZ","QRA","QRB","QRC","QRD","QRE","QRF",
"QRG","QRH","QRI","QRJ","QRK","QRL","QRM","QRN","QRO","QRP","QRQ","QRR","QRS","QRT","QRU","QRV",
"QRW","QRX","QRY","QRZ","QSA","QSB","QSC","QSD","QSE","QSF","QSG","QSH","QSI","QSJ","QSK","QSL",
"QSM","QSN","QSO","QSP","QSQ","QSR","QSS","QST","QSU","QSV","QSW","QSX","QSY","QSZ","QTA","QTB",
"QTC","QTD","QTE","QTF","QTG","QTH","QTI","QTJ","QTK","QTL","QTM","QTN","QTO","QTP","QTQ","QTR",
"QTS","QTT","QTU","QTV","QTW","QTX","QTY","QTZ","QUA","QUB","QUC","QUD","QUE","QUF","QUG","QUH",
"QUI","QUJ","QUK","QUL","QUM","QUN","QUO","QUP","QUQ","QUR","QUS","QUT","QUU","QUV","QUW","QUX",
"QUY","QUZ","QVA","QVB","QVC","QVD","QVE","QVF","QVG","QVH","QVI","QVJ","QVK","QVL","QVM","QVN",
"QVO","QVP","QVQ","QVR","QVS","QVT","QVU","QVV","QVW","QVX","QVY","QVZ","QWA","QWB","QWC","QWD",
"QWE","QWF","QWG","QWH","QWI","QWJ","QWK","QWL","QWM","QWN","QWO","QWP","QWQ","QWR","QWS","QWT",
"QWU","QWV","QWW","QWX","QWY","QWZ","QXA","QXB","QXC","QXD","QXE","QXF","QXG","QXH","QXI","QXJ",
"QXK","QXL","QXM","QXN","QXO","QXP","QXQ","QXR","QXS","QXT","QXU","QXV","QXW","QXX","QXY","QXZ",
"QYA","QYB","QYC","QYD","QYE","QYF","QYG","QYH","QYI","QYJ","QYK","QYL","QYM","QYN","QYO","QYP",
"QYQ","QYR","QYS","QYT","QYU","QYV","QYW","QYX","QYY","QYZ","QZA","QZB","QZC","QZD","QZE","QZF",
"QZG","QZH","QZI","QZJ","QZK","QZL","QZM","QZN","QZO","QZP","QZQ","QZR","QZS","QZT","QZU","QZV",
"QZW","QZX","QZY","QZZ","RAA","RAB","RAC","RAD","RAE","RAF","RAG","RAH","RAI","RAJ","RAK","RAL",
"RAM","RAN","RAO","RAP","RAQ","RAR","RAS","RAT","RAU","RAV","RAW","RAX","RAY","RAZ","RBA","RBB",
"RBC","RBD","RBE","RBF","RBG","RBH","RBI","RBJ","RBK","RBL","RBM","RBN","RBO","RBP","RBQ","RBR",
"RBS","RBT","RBU","RBV","RBW","RBX","RBY","RBZ","RCA","RCB","RCC","RCD","RCE","RCF","RCG","RCH",
"RCI","RCJ","RCK","RCL","RCM","RCN","RCO","RCP","RCQ","RCR","RCS","RCT","RCU","RCV","RCW","RCX",
"RCY","RCZ","RDA","RDB","RDC","RDD","RDE","RDF","RDG","RDH","RDI","RDJ","RDK","RDL","RDM","RDN",
"RDO","RDP","RDQ","RDR","RDS","RDT","RDU","RDV","RDW","RDX","RDY","RDZ","REA","REB","REC","RED",
"REE","REF","REG","REH","REI","REJ","REK","REL","REM","REN","REO","REP","REQ","RER","RES","RET",
"REU","REV","REW","REX","REY","REZ","RFA","RFB","RFC","RFD","RFE","RFF","RFG","RFH","RFI","RFJ",
"RFK","RFL","RFM","RFN","RFO","RFP","RFQ","RFR","RFS","RFT","RFU","RFV","RFW","RFX","RFY","RFZ",
"RGA","RGB","RGC","RGD","RGE","RGF","RGG","RGH","RGI","RGJ","RGK","RGL","RGM","RGN","RGO","RGP",
"RGQ","RGR","RGS","RGT","RGU","RGV","RGW","RGX","RGY","RGZ","RHA","RHB","RHC","RHD","RHE","RHF",
"RHG","RHH","RHI","RHJ","RHK","RHL","RHM","RHN","RHO","RHP","RHQ","RHR","RHS","RHT","RHU","RHV",
"RHW","RHX","RHY","RHZ","RIA","RIB","RIC","RID","RIE","RIF","RIG","RIH","RII","RIJ","RIK","RIL",
"RIM","RIN","RIO","RIP","RIQ","RIR","RIS","RIT","RIU","RIV","RIW","RIX","RIY","RIZ","RJA","RJB",
"RJC","RJD","RJE","RJF","RJG","RJH","RJI","RJJ","RJK","RJL","RJM","RJN","RJO","RJP","RJQ","RJR",
"RJS","RJT","RJU","RJV","RJW","RJX","RJY","RJZ","RKA","RKB","RKC","RKD","RKE","RKF","RKG","RKH",
"RKI","RKJ","RKK","RKL","RKM","RKN","RKO","RKP","RKQ","RKR","RKS","RKT","RKU","RKV","RKW","RKX",
"RKY","RKZ","RLA","RLB","RLC","RLD","RLE","RLF","RLG","RLH","RLI","RLJ","RLK","RLL","RLM","RLN",
"RLO","RLP","RLQ","RLR","RLS","RLT","RLU","RLV","RLW","RLX","RLY","RLZ","RMA","RMB","RMC","RMD",
"RME","RMF","RMG","RMH","RMI","RMJ","RMK","RML","RMM","RMN","RMO","RMP","RMQ","RMR","RMS","RMT",
"RMU","RMV","RMW","RMX","RMY","RMZ","RNA","RNB","RNC","RND","RNE","RNF","RNG","RNH","RNI","RNJ",
"RNK","RNL","RNM","RNN","RNO","RNP","RNQ","RNR","RNS","RNT","RNU","RNV","RNW","RNX","RNY","RNZ",
"ROA","ROB","ROC","ROD","ROE","ROF","ROG","ROH","ROI","ROJ","ROK","ROL","ROM","RON","ROO","ROP",
"ROQ","ROR","ROS","ROT","ROU","ROV","ROW","ROX","ROY","ROZ","RPA","RPB","RPC","RPD","RPE","RPF",
"RPG","RPH","RPI","RPJ","RPK","RPL","RPM","RPN","RPO","RPP","RPQ","RPR","RPS","RPT","RPU","RPV",
"RPW","RPX","RPY","RPZ","RQA","RQB","RQC","RQD","RQE","RQF","RQG","RQH","RQI","RQJ","RQK","RQL",
"RQM","RQN","RQO","RQP","RQQ","RQR","RQS","RQT","RQU","RQV","RQW","RQX","RQY","RQZ","RRA","RRB",
"RRC","RRD","RRE","RRF","RRG","RRH","RRI","RRJ","RRK","RRL","RRM","RRN","RRO","RRP","RRQ","RRR",
"RRS","RRT","RRU","RRV","RRW","RRX","RRY","RRZ","RSA","RSB","RSC","RSD","RSE","RSF","RSG","RSH",
"RSI","RSJ","RSK","RSL","RSM","RSN","RSO","RSP","RSQ","RSR","RSS","RST","RSU","RSV","RSW","RSX",
"RSY","RSZ","RTA","RTB","RTC","RTD","RTE","RTF","RTG","RTH","RTI","RTJ","RTK","RTL","RTM","RTN",
"RTO","RTP","RTQ","RTR","RTS","RTT","RTU","RTV","RTW","RTX","RTY","RTZ","RUA","RUB","RUC","RUD",
"RUE","RUF","RUG","RUH","RUI","RUJ","RUK","RUL","RUM","RUN","RUO","RUP","RUQ","RUR","RUS","RUT",
"RUU","RUV","RUW","RUX","RUY","RUZ","RVA","RVB","RVC","RVD","RVE","RVF","RVG","RVH","RVI","RVJ",
"RVK","RVL","RVM","RVN","RVO","RVP","RVQ","RVR","RVS","RVT","RVU","RVV","RVW","RVX","RVY","RVZ",
"RWA","RWB","RWC","RWD","RWE","RWF","RWG","RWH","RWI","RWJ","RWK","RWL","RWM","RWN","RWO","RWP",
"RWQ","RWR","RWS","RWT","RWU","RWV","RWW","RWX","RWY","RWZ","RXA","RXB","RXC","RXD","RXE","RXF",
"RXG","RXH","RXI","RXJ","RXK","RXL","RXM","RXN","RXO","RXP","RXQ","RXR","RXS","RXT","RXU","RXV",
"RXW","RXX","RXY","RXZ","RYA","RYB","RYC","RYD","RYE","RYF","RYG","RYH","RYI","RYJ","RYK","RYL",
"RYM","RYN","RYO","RYP","RYQ","RYR","RYS","RYT","RYU","RYV","RYW","RYX","RYY","RYZ","RZA","RZB",
"RZC","RZD","RZE","RZF","RZG","RZH","RZI","RZJ","RZK","RZL","RZM","RZN","RZO","RZP","RZQ","RZR",
"RZS","RZT","RZU","RZV","RZW","RZX","RZY","RZZ","SAA","SAB","SAC","SAD","SAE","SAF","SAG","SAH",
"SAI","SAJ","SAK","SAL","SAM","SAN","SAO","SAP","SAQ","SAR","SAS","SAT","SAU","SAV","SAW","SAX",
"SAY","SAZ","SBA","SBB","SBC","SBD","SBE","SBF","SBG","SBH","SBI","SBJ","SBK","SBL","SBM","SBN",
"SBO","SBP","SBQ","SBR","SBS","SBT","SBU","SBV","SBW","SBX","SBY","SBZ","SCA","SCB","SCC","SCD",
"SCE","SCF","SCG","SCH","SCI","SCJ","SCK","SCL","SCM","SCN","SCO","SCP","SCQ","SCR","SCS","SCT",
"SCU","SCV","SCW","SCX","SCY","SCZ","SDA","SDB","SDC","SDD","SDE","SDF","SDG","SDH","SDI","SDJ",
"SDK","SDL","SDM","SDN","SDO","SDP","SDQ","SDR","SDS","SDT","SDU","SDV","SDW","SDX","SDY","SDZ",
"SEA","SEB","SEC","SED","SEE","SEF","SEG","SEH","SEI","SEJ","SEK","SEL","SEM","SEN","SEO","SEP",
"SEQ","SER","SES","SET","SEU","SEV","SEW","SEX","SEY","SEZ","SFA","SFB","SFC","SFD","SFE","SFF",
"SFG","SFH","SFI","SFJ","SFK","SFL","SFM","SFN","SFO","SFP","SFQ","SFR","SFS","SFT","SFU","SFV",
"SFW","SFX","SFY","SFZ","SGA","SGB","SGC","SGD","SGE","SGF","SGG","SGH","SGI","SGJ","SGK","SGL",
"SGM","SGN","SGO","SGP","SGQ","SGR","SGS","SGT","SGU","SGV","SGW","SGX","SGY","SGZ","SHA","SHB",
"SHC","SHD","SHE","SHF","SHG","SHH","SHI","SHJ","SHK","SHL","SHM","SHN","SHO","SHP","SHQ","SHR",
"SHS","SHT","SHU","SHV","SHW","SHX","SHY","SHZ","SIA","SIB","SIC","SID","SIE","SIF","SIG","SIH",
"SII","SIJ","SIK","SIL","SIM","SIN","SIO","SIP","SIQ","SIR","SIS","SIT","SIU","SIV","SIW","SIX",
"SIY","SIZ","SJA","SJB","SJC","SJD","SJE","SJF","SJG","SJH","SJI","SJJ","SJK","SJL","SJM","SJN",
"SJO","SJP","SJQ","SJR","SJS","SJT","SJU","SJV","SJW","SJX","SJY","SJZ","SKA","SKB","SKC","SKD",
"SKE","SKF","SKG","SKH","SKI","SKJ","SKK","SKL","SKM","SKN","SKO","SKP","SKQ","SKR","SKS","SKT",
"SKU","SKV","SKW","SKX","SKY","SKZ","SLA","SLB","SLC","SLD","SLE","SLF","SLG","SLH","SLI","SLJ",
"SLK","SLL","SLM","SLN","SLO","SLP","SLQ","SLR","SLS","SLT","SLU","SLV","SLW","SLX","SLY","SLZ",
"SMA","SMB","SMC","SMD","SME","SMF","SMG","SMH","SMI","SMJ","SMK","SML","SMM","SMN","SMO","SMP",
"SMQ","SMR","SMS","SMT","SMU","SMV","SMW","SMX","SMY","SMZ","SNA","SNB","SNC","SND","SNE","SNF",
"SNG","SNH","SNI","SNJ","SNK","SNL","SNM","SNN","SNO","SNP","SNQ","SNR","SNS","SNT","SNU","SNV",
"SNW","SNX","SNY","SNZ","SOA","SOB","SOC","SOD","SOE","SOF","SOG","SOH","SOI","SOJ","SOK","SOL",
"SOM","SON","SOO","SOP","SOQ","SOR","SOS","SOT","SOU","SOV","SOW","SOX","SOY","SOZ","SPA","SPB",
"SPC","SPD","SPE","SPF","SPG","SPH","SPI","SPJ","SPK","SPL","SPM","SPN","SPO","SPP","SPQ","SPR",
"SPS","SPT","SPU","SPV","SPW","SPX","SPY","SPZ","SQA","SQB","SQC","SQD","SQE","SQF","SQG","SQH",
"SQI","SQJ","SQK","SQL","SQM","SQN","SQO","SQP","SQQ","SQR","SQS","SQT","SQU","SQV","SQW","SQX",
"SQY","SQZ","SRA","SRB","SRC","SRD","SRE","SRF","SRG","SRH","SRI","SRJ","SRK","SRL","SRM","SRN",
"SRO","SRP","SRQ","SRR","SRS","SRT","SRU","SRV","SRW","SRX","SRY","SRZ","SSA","SSB","SSC","SSD",
"SSE","SSF","SSG","SSH","SSI","SSJ","SSK","SSL","SSM","SSN","SSO","SSP","SSQ","SSR","SSS","SST",
"SSU","SSV","SSW","SSX","SSY","SSZ","STA","STB","STC","STD","STE","STF","STG","STH","STI","STJ",
"STK","STL","STM","STN","STO","STP","STQ","STR","STS","STT","STU","STV","STW","STX","STY","STZ",
"SUA","SUB","SUC","SUD","SUE","SUF","SUG","SUH","SUI","SUJ","SUK","SUL","SUM","SUN","SUO","SUP",
"SUQ","SUR","SUS","SUT","SUU","SUV","SUW","SUX","SUY","SUZ","SVA","SVB","SVC","SVD","SVE","SVF",
"SVG","SVH","SVI","SVJ","SVK","SVL","SVM","SVN","SVO","SVP","SVQ","SVR","SVS","SVT","SVU","SVV",
"SVW","SVX","SVY","SVZ","SWA","SWB","SWC","SWD","SWE","SWF","SWG","SWH","SWI","SWJ","SWK","SWL",
"SWM","SWN","SWO","SWP","SWQ","SWR","SWS","SWT","SWU","SWV","SWW","SWX","SWY","SWZ","SXA","SXB",
"SXC","SXD","SXE","SXF","SXG","SXH","SXI","SXJ","SXK","SXL","SXM","SXN","SXO","SXP","SXQ","SXR",
"SXS","SXT","SXU","SXV","SXW","SXX","SXY","SXZ","SYA","SYB","SYC","SYD","SYE","SYF","SYG","SYH",
"SYI","SYJ","SYK","SYL","SYM","SYN","SYO","SYP","SYQ","SYR","SYS","SYT","SYU","SYV","SYW","SYX",
"SYY","SYZ","SZA","SZB","SZC","SZD","SZE","SZF","SZG","SZH","SZI","SZJ","SZK","SZL","SZM","SZN",
"SZO","SZP","SZQ","SZR","SZS","SZT","SZU","SZV","SZW","SZX","SZY","SZZ",
/*
    TAA to TTV - 516 triplets - intentionally omitted
*/
"TTW","TTX","TTY","TTZ","TUA","TUB","TUC","TUD","TUE","TUF","TUG","TUH","TUI","TUJ","TUK","TUL",
"TUM","TUN","TUO","TUP","TUQ","TUR","TUS","TUT","TUU","TUV","TUW","TUX","TUY","TUZ","TVA","TVB",
"TVC","TVD","TVE","TVF","TVG","TVH","TVI","TVJ","TVK","TVL","TVM","TVN","TVO","TVP","TVQ","TVR",
"TVS","TVT","TVU","TVV","TVW","TVX","TVY","TVZ","TWA","TWB","TWC","TWD","TWE","TWF","TWG","TWH",
"TWI","TWJ","TWK","TWL","TWM","TWN","TWO","TWP","TWQ","TWR","TWS","TWT","TWU","TWV","TWW","TWX",
"TWY","TWZ","TXA","TXB","TXC","TXD","TXE","TXF","TXG","TXH","TXI","TXJ","TXK","TXL","TXM","TXN",
"TXO","TXP","TXQ","TXR","TXS","TXT","TXU","TXV","TXW","TXX","TXY","TXZ","TYA","TYB","TYC","TYD",
"TYE","TYF","TYG","TYH","TYI","TYJ","TYK","TYL","TYM","TYN","TYO","TYP","TYQ","TYR","TYS","TYT",
"TYU","TYV","TYW","TYX","TYY","TYZ","TZA","TZB","TZC","TZD","TZE","TZF","TZG","TZH","TZI","TZJ",
"TZK","TZL","TZM","TZN","TZO","TZP","TZQ","TZR","TZS","TZT","TZU","TZV","TZW","TZX","TZY","TZZ",
"UAA","UAB","UAC","UAD","UAE","UAF","UAG","UAH","UAI","UAJ","UAK","UAL","UAM","UAN","UAO","UAP",
"UAQ","UAR","UAS","UAT","UAU","UAV","UAW","UAX","UAY","UAZ","UBA","UBB","UBC","UBD","UBE","UBF",
"UBG","UBH","UBI","UBJ","UBK","UBL","UBM","UBN","UBO","UBP","UBQ","UBR","UBS","UBT","UBU","UBV",
"UBW","UBX","UBY","UBZ","UCA","UCB","UCC","UCD","UCE","UCF","UCG","UCH","UCI","UCJ","UCK","UCL",
"UCM","UCN","UCO","UCP","UCQ","UCR","UCS","UCT","UCU","UCV","UCW","UCX","UCY","UCZ","UDA","UDB",
"UDC","UDD","UDE","UDF","UDG","UDH","UDI","UDJ","UDK","UDL","UDM","UDN","UDO","UDP","UDQ","UDR",
"UDS","UDT","UDU","UDV","UDW","UDX","UDY","UDZ","UEA","UEB","UEC","UED","UEE","UEF","UEG","UEH",
"UEI","UEJ","UEK","UEL","UEM","UEN","UEO","UEP","UEQ","UER","UES","UET","UEU","UEV","UEW","UEX",
"UEY","UEZ","UFA","UFB","UFC","UFD","UFE","UFF","UFG","UFH","UFI","UFJ","UFK","UFL","UFM","UFN",
"UFO","UFP","UFQ","UFR","UFS","UFT","UFU","UFV","UFW","UFX","UFY","UFZ","UGA","UGB","UGC","UGD",
"UGE","UGF","UGG","UGH","UGI","UGJ","UGK","UGL","UGM","UGN","UGO","UGP","UGQ","UGR","UGS","UGT",
"UGU","UGV","UGW","UGX","UGY","UGZ","UHA","UHB","UHC","UHD","UHE","UHF","UHG","UHH","UHI","UHJ",
"UHK","UHL","UHM","UHN","UHO","UHP","UHQ","UHR","UHS","UHT","UHU","UHV","UHW","UHX","UHY","UHZ",
"UIA","UIB","UIC","UID","UIE","UIF","UIG","UIH","UII","UIJ","UIK","UIL","UIM","UIN","UIO","UIP",
"UIQ","UIR","UIS","UIT","UIU","UIV","UIW","UIX","UIY","UIZ","UJA","UJB","UJC","UJD","UJE","UJF",
"UJG","UJH","UJI","UJJ","UJK","UJL","UJM","UJN","UJO","UJP","UJQ","UJR","UJS","UJT","UJU","UJV",
"UJW","UJX","UJY","UJZ","UKA","UKB","UKC","UKD","UKE","UKF","UKG","UKH","UKI","UKJ","UKK","UKL",
"UKM","UKN","UKO","UKP","UKQ","UKR","UKS","UKT","UKU","UKV","UKW","UKX","UKY","UKZ","ULA","ULB",
"ULC","ULD","ULE","ULF","ULG","ULH","ULI","ULJ","ULK","ULL","ULM","ULN","ULO","ULP","ULQ","ULR",
"ULS","ULT","ULU","ULV","ULW","ULX","ULY","ULZ","UMA","UMB","UMC","UMD","UME","UMF","UMG","UMH",
"UMI","UMJ","UMK","UML","UMM","UMN","UMO","UMP","UMQ","UMR","UMS","UMT","UMU","UMV","UMW","UMX",
"UMY","UMZ","UNA","UNB","UNC","UND","UNE","UNF","UNG","UNH","UNI","UNJ","UNK","UNL","UNM","UNN",
"UNO","UNP","UNQ","UNR","UNS","UNT","UNU","UNV","UNW","UNX","UNY","UNZ","UOA","UOB","UOC","UOD",
"UOE","UOF","UOG","UOH","UOI","UOJ","UOK","UOL","UOM","UON","UOO","UOP","UOQ","UOR","UOS","UOT",
"UOU","UOV","UOW","UOX","UOY","UOZ","UPA","UPB","UPC","UPD","UPE","UPF","UPG","UPH","UPI","UPJ",
"UPK","UPL","UPM","UPN","UPO","UPP","UPQ","UPR","UPS","UPT","UPU","UPV","UPW","UPX","UPY","UPZ",
"UQA","UQB","UQC","UQD","UQE","UQF","UQG","UQH","UQI","UQJ","UQK","UQL","UQM","UQN","UQO","UQP",
"UQQ","UQR","UQS","UQT","UQU","UQV","UQW","UQX","UQY","UQZ","URA","URB","URC","URD","URE","URF",
"URG","URH","URI","URJ","URK","URL","URM","URN","URO","URP","URQ","URR","URS","URT","URU","URV",
"URW","URX","URY","URZ","USA","USB","USC","USD","USE","USF","USG","USH","USI","USJ","USK","USL",
"USM","USN","USO","USP","USQ","USR","USS","UST","USU","USV","USW","USX","USY","USZ","UTA","UTB",
"UTC","UTD","UTE","UTF","UTG","UTH","UTI","UTJ","UTK","UTL","UTM","UTN","UTO","UTP","UTQ","UTR",
"UTS","UTT","UTU","UTV","UTW","UTX","UTY","UTZ","UUA","UUB","UUC","UUD","UUE","UUF","UUG","UUH",
"UUI","UUJ","UUK","UUL","UUM","UUN","UUO","UUP","UUQ","UUR","UUS","UUT","UUU","UUV","UUW","UUX",
"UUY","UUZ","UVA","UVB","UVC","UVD","UVE","UVF","UVG","UVH","UVI","UVJ","UVK","UVL","UVM","UVN",
"UVO","UVP","UVQ","UVR","UVS","UVT","UVU","UVV","UVW","UVX","UVY","UVZ","UWA","UWB","UWC","UWD",
"UWE","UWF","UWG","UWH","UWI","UWJ","UWK","UWL","UWM","UWN","UWO","UWP","UWQ","UWR","UWS","UWT",
"UWU","UWV","UWW","UWX","UWY","UWZ","UXA","UXB","UXC","UXD","UXE","UXF","UXG","UXH","UXI","UXJ",
"UXK","UXL","UXM","UXN","UXO","UXP","UXQ","UXR","UXS","UXT","UXU","UXV","UXW","UXX","UXY","UXZ",
"UYA","UYB","UYC","UYD","UYE","UYF","UYG","UYH","UYI","UYJ","UYK","UYL","UYM","UYN","UYO","UYP",
"UYQ","UYR","UYS","UYT","UYU","UYV","UYW","UYX","UYY","UYZ","UZA","UZB","UZC","UZD","UZE","UZF",
"UZG","UZH","UZI","UZJ","UZK","UZL","UZM","UZN","UZO","UZP","UZQ","UZR","UZS","UZT","UZU","UZV",
"UZW","UZX","UZY","UZZ","VAA","VAB","VAC","VAD","VAE","VAF","VAG","VAH","VAI","VAJ","VAK","VAL",
"VAM","VAN","VAO","VAP","VAQ","VAR","VAS","VAT","VAU","VAV","VAW","VAX","VAY","VAZ","VBA","VBB",
"VBC","VBD","VBE","VBF","VBG","VBH","VBI","VBJ","VBK","VBL","VBM","VBN","VBO","VBP","VBQ","VBR",
"VBS","VBT","VBU","VBV","VBW","VBX","VBY","VBZ","VCA","VCB","VCC","VCD","VCE","VCF","VCG","VCH",
"VCI","VCJ","VCK","VCL","VCM","VCN","VCO","VCP","VCQ","VCR","VCS","VCT","VCU","VCV","VCW","VCX",
"VCY","VCZ","VDA","VDB","VDC","VDD","VDE","VDF","VDG","VDH","VDI","VDJ","VDK","VDL","VDM","VDN",
"VDO","VDP","VDQ","VDR","VDS","VDT","VDU","VDV","VDW","VDX","VDY","VDZ","VEA","VEB","VEC","VED",
"VEE","VEF","VEG","VEH","VEI","VEJ","VEK","VEL","VEM","VEN","VEO","VEP","VEQ","VER","VES","VET",
"VEU","VEV","VEW","VEX","VEY","VEZ","VFA","VFB","VFC","VFD","VFE","VFF","VFG","VFH","VFI","VFJ",
"VFK","VFL","VFM","VFN","VFO","VFP","VFQ","VFR","VFS","VFT","VFU","VFV","VFW","VFX","VFY","VFZ",
"VGA","VGB","VGC","VGD","VGE","VGF","VGG","VGH","VGI","VGJ","VGK","VGL","VGM","VGN","VGO","VGP",
"VGQ","VGR","VGS","VGT","VGU","VGV","VGW","VGX","VGY","VGZ","VHA","VHB","VHC","VHD","VHE","VHF",
"VHG","VHH","VHI","VHJ","VHK","VHL","VHM","VHN","VHO","VHP","VHQ","VHR","VHS","VHT","VHU","VHV",
"VHW","VHX","VHY","VHZ","VIA","VIB","VIC","VID","VIE","VIF","VIG","VIH","VII","VIJ","VIK","VIL",
"VIM","VIN","VIO","VIP","VIQ","VIR","VIS","VIT","VIU","VIV","VIW","VIX","VIY","VIZ","VJA","VJB",
"VJC","VJD","VJE","VJF","VJG","VJH","VJI","VJJ","VJK","VJL","VJM","VJN","VJO","VJP","VJQ","VJR",
"VJS","VJT","VJU","VJV","VJW","VJX","VJY","VJZ","VKA","VKB","VKC","VKD","VKE","VKF","VKG","VKH",
"VKI","VKJ","VKK","VKL","VKM","VKN","VKO","VKP","VKQ","VKR","VKS","VKT","VKU","VKV","VKW","VKX",
"VKY","VKZ","VLA","VLB","VLC","VLD","VLE","VLF","VLG","VLH","VLI","VLJ","VLK","VLL","VLM","VLN",
"VLO","VLP","VLQ","VLR","VLS","VLT","VLU","VLV","VLW","VLX","VLY","VLZ","VMA","VMB","VMC","VMD",
"VME","VMF","VMG","VMH","VMI","VMJ","VMK","VML","VMM","VMN","VMO","VMP","VMQ","VMR","VMS","VMT",
"VMU","VMV","VMW","VMX","VMY","VMZ","VNA","VNB","VNC","VND","VNE","VNF","VNG","VNH","VNI","VNJ",
"VNK","VNL","VNM","VNN","VNO","VNP","VNQ","VNR","VNS","VNT","VNU","VNV","VNW","VNX","VNY","VNZ",
"VOA","VOB","VOC","VOD","VOE","VOF","VOG","VOH","VOI","VOJ","VOK","VOL","VOM","VON","VOO","VOP",
"VOQ","VOR","VOS","VOT","VOU","VOV","VOW","VOX","VOY","VOZ","VPA","VPB","VPC","VPD","VPE","VPF",
"VPG","VPH","VPI","VPJ","VPK","VPL","VPM","VPN","VPO","VPP","VPQ","VPR","VPS","VPT","VPU","VPV",
"VPW","VPX","VPY","VPZ","VQA","VQB","VQC","VQD","VQE","VQF","VQG","VQH","VQI","VQJ","VQK","VQL",
"VQM","VQN","VQO","VQP","VQQ","VQR","VQS","VQT","VQU","VQV","VQW","VQX","VQY","VQZ","VRA","VRB",
"VRC","VRD","VRE","VRF","VRG","VRH","VRI","VRJ","VRK","VRL","VRM","VRN","VRO","VRP","VRQ","VRR",
"VRS","VRT","VRU","VRV","VRW","VRX","VRY","VRZ","VSA","VSB","VSC","VSD","VSE","VSF","VSG","VSH",
"VSI","VSJ","VSK","VSL","VSM","VSN","VSO","VSP","VSQ","VSR","VSS","VST","VSU","VSV","VSW","VSX",
"VSY","VSZ","VTA","VTB","VTC","VTD","VTE","VTF","VTG","VTH","VTI","VTJ","VTK","VTL","VTM","VTN",
"VTO","VTP","VTQ","VTR","VTS","VTT","VTU","VTV","VTW","VTX","VTY","VTZ","VUA","VUB","VUC","VUD",
"VUE","VUF","VUG","VUH","VUI","VUJ","VUK","VUL","VUM","VUN","VUO","VUP","VUQ","VUR","VUS","VUT",
"VUU","VUV","VUW","VUX","VUY","VUZ","VVA","VVB","VVC","VVD","VVE","VVF","VVG","VVH","VVI","VVJ",
"VVK","VVL","VVM","VVN","VVO","VVP","VVQ","VVR","VVS","VVT","VVU","VVV","VVW","VVX","VVY","VVZ",
"VWA","VWB","VWC","VWD","VWE","VWF","VWG","VWH","VWI","VWJ","VWK","VWL","VWM","VWN","VWO","VWP",
"VWQ","VWR","VWS","VWT","VWU","VWV","VWW","VWX","VWY","VWZ","VXA","VXB","VXC","VXD","VXE","VXF",
"VXG","VXH","VXI","VXJ","VXK","VXL","VXM","VXN","VXO","VXP","VXQ","VXR","VXS","VXT","VXU","VXV",
"VXW","VXX","VXY","VXZ","VYA","VYB","VYC","VYD","VYE","VYF","VYG","VYH","VYI","VYJ","VYK","VYL",
"VYM","VYN","VYO","VYP","VYQ","VYR","VYS","VYT","VYU","VYV","VYW","VYX","VYY","VYZ","VZA","VZB",
"VZC","VZD","VZE","VZF","VZG","VZH","VZI","VZJ","VZK","VZL","VZM","VZN","VZO","VZP","VZQ","VZR",
"VZS","VZT","VZU","VZV","VZW","VZX","VZY","VZZ","WAA","WAB","WAC","WAD","WAE","WAF","WAG","WAH",
"WAI","WAJ","WAK","WAL","WAM","WAN","WAO","WAP","WAQ","WAR","WAS","WAT","WAU","WAV","WAW","WAX",
"WAY","WAZ","WBA","WBB","WBC","WBD","WBE","WBF","WBG","WBH","WBI","WBJ","WBK","WBL","WBM","WBN",
"WBO","WBP","WBQ","WBR","WBS","WBT","WBU","WBV","WBW","WBX","WBY","WBZ","WCA","WCB","WCC","WCD",
"WCE","WCF","WCG","WCH","WCI","WCJ","WCK","WCL","WCM","WCN","WCO","WCP","WCQ","WCR","WCS","WCT",
"WCU","WCV","WCW","WCX","WCY","WCZ","WDA","WDB","WDC","WDD","WDE","WDF","WDG","WDH","WDI","WDJ",
"WDK","WDL","WDM","WDN","WDO","WDP","WDQ","WDR","WDS","WDT","WDU","WDV","WDW","WDX","WDY","WDZ",
"WEA","WEB","WEC","WED","WEE","WEF","WEG","WEH","WEI","WEJ","WEK","WEL","WEM","WEN","WEO","WEP",
"WEQ","WER","WES","WET","WEU","WEV","WEW","WEX","WEY","WEZ","WFA","WFB","WFC","WFD","WFE","WFF",
"WFG","WFH","WFI","WFJ","WFK","WFL","WFM","WFN","WFO","WFP","WFQ","WFR","WFS","WFT","WFU","WFV",
"WFW","WFX","WFY","WFZ","WGA","WGB","WGC","WGD","WGE","WGF","WGG","WGH","WGI","WGJ","WGK","WGL",
"WGM","WGN","WGO","WGP","WGQ","WGR","WGS","WGT","WGU","WGV","WGW","WGX","WGY","WGZ","WHA","WHB",
"WHC","WHD","WHE","WHF","WHG","WHH","WHI","WHJ","WHK","WHL","WHM","WHN","WHO","WHP","WHQ","WHR",
"WHS","WHT","WHU","WHV","WHW","WHX","WHY","WHZ","WIA","WIB","WIC","WID","WIE","WIF","WIG","WIH",
"WII","WIJ","WIK","WIL","WIM","WIN","WIO","WIP","WIQ","WIR","WIS","WIT","WIU","WIV","WIW","WIX",
"WIY","WIZ","WJA","WJB","WJC","WJD","WJE","WJF","WJG","WJH","WJI","WJJ","WJK","WJL","WJM","WJN",
"WJO","WJP","WJQ","WJR","WJS","WJT","WJU","WJV","WJW","WJX","WJY","WJZ","WKA","WKB","WKC","WKD",
"WKE","WKF","WKG","WKH","WKI","WKJ","WKK","WKL","WKM","WKN","WKO","WKP","WKQ","WKR","WKS","WKT",
"WKU","WKV","WKW","WKX","WKY","WKZ","WLA","WLB","WLC","WLD","WLE","WLF","WLG","WLH","WLI","WLJ",
"WLK","WLL","WLM","WLN","WLO","WLP","WLQ","WLR","WLS","WLT","WLU","WLV","WLW","WLX","WLY","WLZ",
"WMA","WMB","WMC","WMD","WME","WMF","WMG","WMH","WMI","WMJ","WMK","WML","WMM","WMN","WMO","WMP",
"WMQ","WMR","WMS","WMT","WMU","WMV","WMW","WMX","WMY","WMZ","WNA","WNB","WNC","WND","WNE","WNF",
"WNG","WNH","WNI","WNJ","WNK","WNL","WNM","WNN","WNO","WNP","WNQ","WNR","WNS","WNT","WNU","WNV",
"WNW","WNX","WNY","WNZ","WOA","WOB","WOC","WOD","WOE","WOF","WOG","WOH","WOI","WOJ","WOK","WOL",
"WOM","WON","WOO","WOP","WOQ","WOR","WOS","WOT","WOU","WOV","WOW","WOX","WOY","WOZ","WPA","WPB",
"WPC","WPD","WPE","WPF","WPG","WPH","WPI","WPJ","WPK","WPL","WPM","WPN","WPO","WPP","WPQ","WPR",
"WPS","WPT","WPU","WPV","WPW","WPX","WPY","WPZ","WQA","WQB","WQC","WQD","WQE","WQF","WQG","WQH",
"WQI","WQJ","WQK","WQL","WQM","WQN","WQO","WQP","WQQ","WQR","WQS","WQT","WQU","WQV","WQW","WQX",
"WQY","WQZ","WRA","WRB","WRC","WRD","WRE","WRF","WRG","WRH","WRI","WRJ","WRK","WRL","WRM","WRN",
"WRO","WRP","WRQ","WRR","WRS","WRT","WRU","WRV","WRW","WRX","WRY","WRZ","WSA","WSB","WSC","WSD",
"WSE","WSF","WSG","WSH","WSI","WSJ","WSK","WSL","WSM","WSN","WSO","WSP","WSQ","WSR","WSS","WST",
"WSU","WSV","WSW","WSX","WSY","WSZ","WTA","WTB","WTC","WTD","WTE","WTF","WTG","WTH","WTI","WTJ",
"WTK","WTL","WTM","WTN","WTO","WTP","WTQ","WTR","WTS","WTT","WTU","WTV","WTW","WTX","WTY","WTZ",
"WUA","WUB","WUC","WUD","WUE","WUF","WUG","WUH","WUI","WUJ","WUK","WUL","WUM","WUN","WUO","WUP",
"WUQ","WUR","WUS","WUT","WUU","WUV","WUW","WUX","WUY","WUZ","WVA","WVB","WVC","WVD","WVE","WVF",
"WVG","WVH","WVI","WVJ","WVK","WVL","WVM","WVN","WVO","WVP","WVQ","WVR","WVS","WVT","WVU","WVV",
"WVW","WVX","WVY","WVZ","WWA","WWB","WWC","WWD","WWE","WWF","WWG","WWH","WWI","WWJ","WWK","WWL",
"WWM","WWN","WWO","WWP","WWQ","WWR","WWS","WWT","WWU","WWV","WWW","WWX","WWY","WWZ","WXA","WXB",
"WXC","WXD","WXE","WXF","WXG","WXH","WXI","WXJ","WXK","WXL","WXM","WXN","WXO","WXP","WXQ","WXR",
"WXS","WXT","WXU","WXV","WXW","WXX","WXY","WXZ","WYA","WYB","WYC","WYD","WYE","WYF","WYG","WYH",
"WYI","WYJ","WYK","WYL","WYM","WYN","WYO","WYP","WYQ","WYR","WYS","WYT","WYU","WYV","WYW","WYX",
"WYY","WYZ","WZA","WZB","WZC","WZD","WZE","WZF","WZG","WZH","WZI","WZJ","WZK","WZL","WZM","WZN",
"WZO","WZP","WZQ","WZR","WZS","WZT","WZU","WZV","WZW","WZX","WZY","WZZ","XAA","XAB","XAC","XAD",
"XAE","XAF","XAG","XAH","XAI","XAJ","XAK","XAL","XAM","XAN","XAO","XAP","XAQ","XAR","XAS","XAT",
"XAU","XAV","XAW","XAX","XAY","XAZ","XBA","XBB","XBC","XBD","XBE","XBF","XBG","XBH","XBI","XBJ",
"XBK","XBL","XBM","XBN","XBO","XBP","XBQ","XBR","XBS","XBT","XBU","XBV","XBW","XBX","XBY","XBZ",
"XCA","XCB","XCC","XCD","XCE","XCF","XCG","XCH","XCI","XCJ","XCK","XCL","XCM","XCN","XCO","XCP",
"XCQ","XCR","XCS","XCT","XCU","XCV","XCW","XCX","XCY","XCZ","XDA","XDB","XDC","XDD","XDE","XDF",
"XDG","XDH","XDI","XDJ","XDK","XDL","XDM","XDN","XDO","XDP","XDQ","XDR","XDS","XDT","XDU","XDV",
"XDW","XDX","XDY","XDZ","XEA","XEB","XEC","XED","XEE","XEF","XEG","XEH","XEI","XEJ","XEK","XEL",
"XEM","XEN","XEO","XEP","XEQ","XER","XES","XET","XEU","XEV","XEW","XEX","XEY","XEZ","XFA","XFB",
"XFC","XFD","XFE","XFF","XFG","XFH","XFI","XFJ","XFK","XFL","XFM","XFN","XFO","XFP","XFQ","XFR",
"XFS","XFT","XFU","XFV","XFW","XFX","XFY","XFZ","XGA","XGB","XGC","XGD","XGE","XGF","XGG","XGH",
"XGI","XGJ","XGK","XGL","XGM","XGN","XGO","XGP","XGQ","XGR","XGS","XGT","XGU","XGV","XGW","XGX",
"XGY","XGZ","XHA","XHB","XHC","XHD","XHE","XHF","XHG","XHH","XHI","XHJ","XHK","XHL","XHM","XHN",
"XHO","XHP","XHQ","XHR","XHS","XHT","XHU","XHV","XHW","XHX","XHY","XHZ","XIA","XIB","XIC","XID",
"XIE","XIF","XIG","XIH","XII","XIJ","XIK","XIL","XIM","XIN","XIO","XIP","XIQ","XIR","XIS","XIT",
"XIU","XIV","XIW","XIX","XIY","XIZ","XJA","XJB","XJC","XJD","XJE","XJF","XJG","XJH","XJI","XJJ",
"XJK","XJL","XJM","XJN","XJO","XJP","XJQ","XJR","XJS","XJT","XJU","XJV","XJW","XJX","XJY","XJZ",
"XKA","XKB","XKC","XKD","XKE","XKF","XKG","XKH","XKI","XKJ","XKK","XKL","XKM","XKN","XKO","XKP",
"XKQ","XKR","XKS","XKT","XKU","XKV","XKW","XKX","XKY","XKZ","XLA","XLB","XLC","XLD","XLE","XLF",
"XLG","XLH","XLI","XLJ","XLK","XLL","XLM","XLN","XLO","XLP","XLQ","XLR","XLS","XLT","XLU","XLV",
"XLW","XLX","XLY","XLZ","XMA","XMB","XMC","XMD","XME","XMF","XMG","XMH","XMI","XMJ","XMK","XML",
"XMM","XMN","XMO","XMP","XMQ","XMR","XMS","XMT","XMU","XMV","XMW","XMX","XMY","XMZ","XNA","XNB",
"XNC","XND","XNE","XNF","XNG","XNH","XNI","XNJ","XNK","XNL","XNM","XNN","XNO","XNP","XNQ","XNR",
"XNS","XNT","XNU","XNV","XNW","XNX","XNY","XNZ","XOA","XOB","XOC","XOD","XOE","XOF","XOG","XOH",
"XOI","XOJ","XOK","XOL","XOM","XON","XOO","XOP","XOQ","XOR","XOS","XOT","XOU","XOV","XOW","XOX",
"XOY","XOZ","XPA","XPB","XPC","XPD","XPE","XPF","XPG","XPH","XPI","XPJ","XPK","XPL","XPM","XPN",
"XPO","XPP","XPQ","XPR","XPS","XPT","XPU","XPV","XPW","XPX","XPY","XPZ","XQA","XQB","XQC","XQD",
"XQE","XQF","XQG","XQH","XQI","XQJ","XQK","XQL","XQM","XQN","XQO","XQP","XQQ","XQR","XQS","XQT",
"XQU","XQV","XQW","XQX","XQY","XQZ","XRA","XRB","XRC","XRD","XRE","XRF","XRG","XRH","XRI","XRJ",
"XRK","XRL","XRM","XRN","XRO","XRP","XRQ","XRR","XRS","XRT","XRU","XRV","XRW","XRX","XRY","XRZ",
"XSA","XSB","XSC","XSD","XSE","XSF","XSG","XSH","XSI","XSJ","XSK","XSL","XSM","XSN","XSO","XSP",
"XSQ","XSR","XSS","XST","XSU","XSV","XSW","XSX","XSY","XSZ","XTA","XTB","XTC","XTD","XTE","XTF",
"XTG","XTH","XTI","XTJ","XTK","XTL","XTM","XTN","XTO","XTP","XTQ","XTR","XTS","XTT","XTU","XTV",
"XTW","XTX","XTY","XTZ","XUA","XUB","XUC","XUD","XUE","XUF","XUG","XUH","XUI","XUJ","XUK","XUL",
"XUM","XUN","XUO","XUP","XUQ","XUR","XUS","XUT","XUU","XUV","XUW","XUX","XUY","XUZ","XVA","XVB",
"XVC","XVD","XVE","XVF","XVG","XVH","XVI","XVJ","XVK","XVL","XVM","XVN","XVO","XVP","XVQ","XVR",
"XVS","XVT","XVU","XVV","XVW","XVX","XVY","XVZ","XWA","XWB","XWC","XWD","XWE","XWF","XWG","XWH",
"XWI","XWJ","XWK","XWL","XWM","XWN","XWO","XWP","XWQ","XWR","XWS","XWT","XWU","XWV","XWW","XWX",
"XWY","XWZ","XXA","XXB","XXC","XXD","XXE","XXF","XXG","XXH","XXI","XXJ","XXK","XXL","XXM","XXN",
"XXO","XXP","XXQ","XXR","XXS","XXT","XXU","XXV","XXW","XXX","XXY","XXZ","XYA","XYB","XYC","XYD",
"XYE","XYF","XYG","XYH","XYI","XYJ","XYK","XYL","XYM","XYN","XYO","XYP","XYQ","XYR","XYS","XYT",
"XYU","XYV","XYW","XYX","XYY","XYZ","XZA","XZB","XZC","XZD","XZE","XZF","XZG","XZH","XZI","XZJ",
"XZK","XZL","XZM","XZN","XZO","XZP","XZQ","XZR","XZS","XZT","XZU","XZV","XZW","XZX","XZY","XZZ",

"YAA","YAB","YAC","YAD",
"YAE","YAF","YAG","YAH","YAI","YAJ","YAK","YAL","YAM","YAN","YAO","YAP","YAQ","YAR","YAS","YAT",
"YAU","YAV","YAW","YAX","YAY","YAZ","YBA","YBB","YBC","YBD","YBE","YBF","YBG","YBH","YBI","YBJ",
"YBK","YBL","YBM","YBN","YBO","YBP","YBQ","YBR","YBS","YBT","YBU","YBV","YBW","YBX","YBY","YBZ",
"YCA","YCB","YCC","YCD","YCE","YCF","YCG","YCH","YCI","YCJ","YCK","YCL","YCM","YCN","YCO","YCP",
"YCQ","YCR","YCS","YCT","YCU","YCV","YCW","YCX","YCY","YCZ","YDA","YDB","YDC","YDD","YDE","YDF",
"YDG","YDH","YDI","YDJ","YDK","YDL","YDM","YDN","YDO","YDP","YDQ","YDR","YDS","YDT","YDU","YDV",
"YDW","YDX","YDY","YDZ","YEA","YEB","YEC","YED","YEE","YEF","YEG","YEH","YEI","YEJ","YEK","YEL",
"YEM","YEN","YEO","YEP","YEQ","YER","YES","YET","YEU","YEV","YEW","YEX","YEY","YEZ","YFA","YFB",
"YFC","YFD","YFE","YFF","YFG","YFH","YFI","YFJ","YFK","YFL","YFM","YFN","YFO","YFP","YFQ","YFR",
"YFS","YFT","YFU","YFV","YFW","YFX","YFY","YFZ","YGA","YGB","YGC","YGD","YGE","YGF","YGG","YGH",
"YGI","YGJ","YGK","YGL","YGM","YGN","YGO","YGP","YGQ","YGR","YGS","YGT","YGU","YGV","YGW","YGX",
"YGY","YGZ","YHA","YHB","YHC","YHD","YHE","YHF","YHG","YHH","YHI","YHJ","YHK","YHL","YHM","YHN",
"YHO","YHP","YHQ","YHR","YHS","YHT","YHU","YHV","YHW","YHX","YHY","YHZ","YIA","YIB","YIC","YID",
"YIE","YIF","YIG","YIH","YII","YIJ","YIK","YIL","YIM","YIN","YIO","YIP","YIQ","YIR","YIS","YIT",
"YIU","YIV","YIW","YIX","YIY","YIZ","YJA","YJB","YJC","YJD","YJE","YJF","YJG","YJH","YJI","YJJ",
"YJK","YJL","YJM","YJN","YJO","YJP","YJQ","YJR","YJS","YJT","YJU","YJV","YJW","YJX","YJY","YJZ",
"YKA","YKB","YKC","YKD","YKE","YKF","YKG","YKH","YKI","YKJ","YKK","YKL","YKM","YKN","YKO","YKP",
"YKQ","YKR","YKS","YKT","YKU","YKV","YKW","YKX","YKY","YKZ","YLA","YLB","YLC","YLD","YLE","YLF",
"YLG","YLH","YLI","YLJ","YLK","YLL","YLM","YLN","YLO","YLP","YLQ","YLR","YLS","YLT","YLU","YLV",
"YLW","YLX","YLY","YLZ","YMA","YMB","YMC","YMD","YME","YMF","YMG","YMH","YMI","YMJ","YMK","YML",
"YMM","YMN","YMO","YMP","YMQ","YMR","YMS","YMT","YMU","YMV","YMW","YMX","YMY","YMZ","YNA","YNB",
"YNC","YND","YNE","YNF","YNG","YNH","YNI","YNJ","YNK","YNL","YNM","YNN","YNO","YNP","YNQ","YNR",
"YNS","YNT","YNU","YNV","YNW","YNX","YNY","YNZ","YOA","YOB","YOC","YOD","YOE","YOF","YOG","YOH",
"YOI","YOJ","YOK","YOL","YOM","YON","YOO","YOP","YOQ","YOR","YOS","YOT","YOU","YOV","YOW","YOX",
"YOY","YOZ","YPA","YPB","YPC","YPD","YPE","YPF","YPG","YPH","YPI","YPJ","YPK","YPL","YPM","YPN",
"YPO","YPP","YPQ","YPR","YPS","YPT","YPU","YPV","YPW","YPX","YPY","YPZ","YQA","YQB","YQC","YQD",
"YQE","YQF","YQG","YQH","YQI","YQJ","YQK","YQL","YQM","YQN","YQO","YQP","YQQ","YQR","YQS","YQT",
"YQU","YQV","YQW","YQX","YQY","YQZ","YRA","YRB","YRC","YRD","YRE","YRF","YRG","YRH","YRI","YRJ",
"YRK","YRL","YRM","YRN","YRO","YRP","YRQ","YRR","YRS","YRT","YRU","YRV","YRW","YRX","YRY","YRZ",
"YSA","YSB","YSC","YSD","YSE","YSF","YSG","YSH","YSI","YSJ","YSK","YSL","YSM","YSN","YSO","YSP",
"YSQ","YSR","YSS","YST","YSU","YSV","YSW","YSX","YSY","YSZ","YTA","YTB","YTC","YTD","YTE","YTF",
"YTG","YTH","YTI","YTJ","YTK","YTL","YTM","YTN","YTO","YTP","YTQ","YTR","YTS","YTT","YTU","YTV",
"YTW","YTX","YTY","YTZ","YUA","YUB","YUC","YUD","YUE","YUF","YUG","YUH","YUI","YUJ","YUK","YUL",
"YUM","YUN","YUO","YUP","YUQ","YUR","YUS","YUT","YUU","YUV","YUW","YUX","YUY","YUZ","YVA","YVB",
"YVC","YVD","YVE","YVF","YVG","YVH","YVI","YVJ","YVK","YVL","YVM","YVN","YVO","YVP","YVQ","YVR",
"YVS","YVT","YVU","YVV","YVW","YVX","YVY","YVZ","YWA","YWB","YWC","YWD","YWE","YWF","YWG","YWH",
"YWI","YWJ","YWK","YWL","YWM","YWN","YWO","YWP","YWQ","YWR","YWS","YWT","YWU","YWV","YWW","YWX",
"YWY","YWZ","YXA","YXB","YXC","YXD","YXE","YXF","YXG","YXH","YXI","YXJ","YXK","YXL","YXM","YXN",
"YXO","YXP","YXQ","YXR","YXS","YXT","YXU","YXV","YXW","YXX","YXY","YXZ","YYA","YYB","YYC","YYD",
"YYE","YYF","YYG","YYH","YYI","YYJ","YYK","YYL","YYM","YYN","YYO","YYP","YYQ","YYR","YYS","YYT",
"YYU","YYV","YYW","YYX","YYY","YYZ","YZA","YZB","YZC","YZD","YZE","YZF","YZG","YZH","YZI","YZJ",
"YZK","YZL","YZM","YZN","YZO","YZP","YZQ","YZR","YZS","YZT","YZU","YZV","YZW","YZX","YZY","YZZ",

"ZAA","ZAB","ZAC","ZAD",
"ZAE","ZAF","ZAG","ZAH","ZAI","ZAJ","ZAK","ZAL","ZAM","ZAN","ZAO","ZAP","ZAQ","ZAR","ZAS","ZAT",
"ZAU","ZAV","ZAW","ZAX","ZAY","ZAZ","ZBA","ZBB","ZBC","ZBD","ZBE","ZBF","ZBG","ZBH","ZBI","ZBJ",
"ZBK","ZBL","ZBM","ZBN","ZBO","ZBP","ZBQ","ZBR","ZBS","ZBT","ZBU","ZBV","ZBW","ZBX","ZBY","ZBZ",
"ZCA","ZCB","ZCC","ZCD","ZCE","ZCF","ZCG","ZCH","ZCI","ZCJ","ZCK","ZCL","ZCM","ZCN","ZCO","ZCP",
"ZCQ","ZCR","ZCS","ZCT","ZCU","ZCV","ZCW","ZCX","ZCY","ZCZ","ZDA","ZDB","ZDC","ZDD","ZDE","ZDF",
"ZDG","ZDH","ZDI","ZDJ","ZDK","ZDL","ZDM","ZDN","ZDO","ZDP","ZDQ","ZDR","ZDS","ZDT","ZDU","ZDV",
"ZDW","ZDX","ZDY","ZDZ","ZEA","ZEB","ZEC","ZED","ZEE","ZEF","ZEG","ZEH","ZEI","ZEJ","ZEK","ZEL",
"ZEM","ZEN","ZEO","ZEP","ZEQ","ZER","ZES","ZET","ZEU","ZEV","ZEW","ZEX","ZEY","ZEZ","ZFA","ZFB",
"ZFC","ZFD","ZFE","ZFF","ZFG","ZFH","ZFI","ZFJ","ZFK","ZFL","ZFM","ZFN","ZFO","ZFP","ZFQ","ZFR",
"ZFS","ZFT","ZFU","ZFV","ZFW","ZFX","ZFY","ZFZ","ZGA","ZGB","ZGC","ZGD","ZGE","ZGF","ZGG","ZGH",
"ZGI","ZGJ","ZGK","ZGL","ZGM","ZGN","ZGO","ZGP","ZGQ","ZGR","ZGS","ZGT","ZGU","ZGV","ZGW","ZGX",
"ZGY","ZGZ","ZHA","ZHB","ZHC","ZHD","ZHE","ZHF","ZHG","ZHH","ZHI","ZHJ","ZHK","ZHL","ZHM","ZHN",
"ZHO","ZHP","ZHQ","ZHR","ZHS","ZHT","ZHU","ZHV","ZHW","ZHX","ZHY","ZHZ","ZIA","ZIB","ZIC","ZID",
"ZIE","ZIF","ZIG","ZIH","ZII","ZIJ","ZIK","ZIL","ZIM","ZIN","ZIO","ZIP","ZIQ","ZIR","ZIS","ZIT",
"ZIU","ZIV","ZIW","ZIX","ZIY","ZIZ","ZJA","ZJB","ZJC","ZJD","ZJE","ZJF","ZJG","ZJH","ZJI","ZJJ",
"ZJK","ZJL","ZJM","ZJN","ZJO","ZJP","ZJQ","ZJR","ZJS","ZJT","ZJU","ZJV","ZJW","ZJX","ZJY","ZJZ",
"ZKA","ZKB","ZKC","ZKD","ZKE","ZKF","ZKG","ZKH","ZKI","ZKJ","ZKK","ZKL","ZKM","ZKN","ZKO","ZKP",
"ZKQ","ZKR","ZKS","ZKT","ZKU","ZKV","ZKW","ZKX","ZKY","ZKZ","ZLA","ZLB","ZLC","ZLD","ZLE","ZLF",
"ZLG","ZLH","ZLI","ZLJ","ZLK","ZLL","ZLM","ZLN","ZLO","ZLP","ZLQ","ZLR","ZLS","ZLT","ZLU","ZLV",
"ZLW","ZLX","ZLY","ZLZ","ZMA","ZMB","ZMC","ZMD","ZME","ZMF","ZMG","ZMH","ZMI","ZMJ","ZMK","ZML",
"ZMM","ZMN","ZMO","ZMP","ZMQ","ZMR","ZMS","ZMT","ZMU","ZMV","ZMW","ZMX","ZMY","ZMZ","ZNA","ZNB",
"ZNC","ZND","ZNE","ZNF","ZNG","ZNH","ZNI","ZNJ","ZNK","ZNL","ZNM","ZNN","ZNO","ZNP","ZNQ","ZNR",
"ZNS","ZNT","ZNU","ZNV","ZNW","ZNX","ZNY","ZNZ","ZOA","ZOB","ZOC","ZOD","ZOE","ZOF","ZOG","ZOH",
"ZOI","ZOJ","ZOK","ZOL","ZOM","ZON","ZOO","ZOP","ZOQ","ZOR","ZOS","ZOT","ZOU","ZOV","ZOW","ZOX",
"ZOY","ZOZ","ZPA","ZPB","ZPC","ZPD","ZPE","ZPF","ZPG","ZPH","ZPI","ZPJ","ZPK","ZPL","ZPM","ZPN",
"ZPO","ZPP","ZPQ","ZPR","ZPS","ZPT","ZPU","ZPV","ZPW","ZPX","ZPY","ZPZ","ZQA","ZQB","ZQC","ZQD",
"ZQE","ZQF","ZQG","ZQH","ZQI","ZQJ","ZQK","ZQL","ZQM","ZQN","ZQO","ZQP","ZQQ","ZQR","ZQS","ZQT",
"ZQU","ZQV","ZQW","ZQX","ZQY","ZQZ","ZRA","ZRB","ZRC","ZRD","ZRE","ZRF","ZRG","ZRH","ZRI","ZRJ",
"ZRK","ZRL","ZRM","ZRN","ZRO","ZRP","ZRQ","ZRR","ZRS","ZRT","ZRU","ZRV","ZRW","ZRX","ZRY","ZRZ",
"ZSA","ZSB","ZSC","ZSD","ZSE","ZSF","ZSG","ZSH","ZSI","ZSJ","ZSK","ZSL","ZSM","ZSN","ZSO","ZSP",
"ZSQ","ZSR","ZSS","ZST","ZSU","ZSV","ZSW","ZSX","ZSY","ZSZ","ZTA","ZTB","ZTC","ZTD","ZTE","ZTF",
"ZTG","ZTH","ZTI","ZTJ","ZTK","ZTL","ZTM","ZTN","ZTO","ZTP","ZTQ","ZTR","ZTS","ZTT","ZTU","ZTV",
"ZTW","ZTX","ZTY","ZTZ","ZUA","ZUB","ZUC","ZUD","ZUE","ZUF","ZUG","ZUH","ZUI","ZUJ","ZUK","ZUL",
"ZUM","ZUN","ZUO","ZUP","ZUQ","ZUR","ZUS","ZUT","ZUU","ZUV","ZUW","ZUX","ZUY","ZUZ","ZVA","ZVB",
"ZVC","ZVD","ZVE","ZVF","ZVG","ZVH","ZVI","ZVJ","ZVK","ZVL","ZVM","ZVN","ZVO","ZVP","ZVQ","ZVR",
"ZVS","ZVT","ZVU","ZVV","ZVW","ZVX","ZVY","ZVZ","ZWA","ZWB","ZWC","ZWD","ZWE","ZWF","ZWG","ZWH",
"ZWI","ZWJ","ZWK","ZWL","ZWM","ZWN","ZWO","ZWP","ZWQ","ZWR","ZWS","ZWT","ZWU","ZWV","ZWW","ZWX",
"ZWY","ZWZ","ZXA","ZXB","ZXC","ZXD","ZXE","ZXF","ZXG","ZXH","ZXI","ZXJ","ZXK","ZXL","ZXM","ZXN",
"ZXO","ZXP","ZXQ","ZXR","ZXS","ZXT","ZXU","ZXV","ZXW","ZXX","ZXY","ZXZ","ZYA","ZYB","ZYC","ZYD",
"ZYE","ZYF","ZYG","ZYH","ZYI","ZYJ","ZYK","ZYL","ZYM","ZYN","ZYO","ZYP","ZYQ","ZYR","ZYS","ZYT",
"ZYU","ZYV","ZYW","ZYX","ZYY","ZYZ","ZZA","ZZB","ZZC","ZZD","ZZE","ZZF","ZZG","ZZH","ZZI","ZZJ",
"ZZK","ZZL","ZZM","ZZN","ZZO","ZZP","ZZQ","ZZR","ZZS","ZZT","ZZU","ZZV","ZZW","ZZX","ZZY","ZZZ"
};

/*
    Doublets
*/

static const char d26[][3] =
{
"AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP",
"AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD","BE","BF",
"BG","BH","BI","BJ","BK","BL","BM","BN","BO","BP","BQ","BR","BS","BT","BU","BV",
"BW","BX","BY","BZ","CA","CB","CC","CD","CE","CF","CG","CH","CI","CJ","CK","CL",
"CM","CN","CO","CP","CQ","CR","CS","CT","CU","CV","CW","CX","CY","CZ","DA","DB",
"DC","DD","DE","DF","DG","DH","DI","DJ","DK","DL","DM","DN","DO","DP","DQ","DR",
"DS","DT","DU","DV","DW","DX","DY","DZ","EA","EB","EC","ED","EE","EF","EG","EH",
"EI","EJ","EK","EL","EM","EN","EO","EP","EQ","ER","ES","ET","EU","EV","EW","EX",
"EY","EZ","FA","FB","FC","FD","FE","FF","FG","FH","FI","FJ","FK","FL","FM","FN",
"FO","FP","FQ","FR","FS","FT","FU","FV","FW","FX","FY","FZ","GA","GB","GC","GD",
"GE","GF","GG","GH","GI","GJ","GK","GL","GM","GN","GO","GP","GQ","GR","GS","GT",
"GU","GV","GW","GX","GY","GZ","HA","HB","HC","HD","HE","HF","HG","HH","HI","HJ",
"HK","HL","HM","HN","HO","HP","HQ","HR","HS","HT","HU","HV","HW","HX","HY","HZ",
"IA","IB","IC","ID","IE","IF","IG","IH","II","IJ","IK","IL","IM","IN","IO","IP",
"IQ","IR","IS","IT","IU","IV","IW","IX","IY","IZ","JA","JB","JC","JD","JE","JF",
"JG","JH","JI","JJ","JK","JL","JM","JN","JO","JP","JQ","JR","JS","JT","JU","JV",
"JW","JX","JY","JZ","KA","KB","KC","KD","KE","KF","KG","KH","KI","KJ","KK","KL",
"KM","KN","KO","KP","KQ","KR","KS","KT","KU","KV","KW","KX","KY","KZ","LA","LB",
"LC","LD","LE","LF","LG","LH","LI","LJ","LK","LL","LM","LN","LO","LP","LQ","LR",
"LS","LT","LU","LV","LW","LX","LY","LZ","MA","MB","MC","MD","ME","MF","MG","MH",
"MI","MJ","MK","ML","MM","MN","MO","MP","MQ","MR","MS","MT","MU","MV","MW","MX",
"MY","MZ","NA","NB","NC","ND","NE","NF","NG","NH","NI","NJ","NK","NL","NM","NN",
"NO","NP","NQ","NR","NS","NT","NU","NV","NW","NX","NY","NZ","OA","OB","OC","OD",
"OE","OF","OG","OH","OI","OJ","OK","OL","OM","ON","OO","OP","OQ","OR","OS","OT",
"OU","OV","OW","OX","OY","OZ","PA","PB","PC","PD","PE","PF","PG","PH","PI","PJ",
"PK","PL","PM","PN","PO","PP","PQ","PR","PS","PT","PU","PV","PW","PX","PY","PZ",
"QA","QB","QC","QD","QE","QF","QG","QH","QI","QJ","QK","QL","QM","QN","QO","QP",
"QQ","QR","QS","QT","QU","QV","QW","QX","QY","QZ","RA","RB","RC","RD","RE","RF",
"RG","RH","RI","RJ","RK","RL","RM","RN","RO","RP","RQ","RR","RS","RT","RU","RV",
"RW","RX","RY","RZ","SA","SB","SC","SD","SE","SF","SG","SH","SI","SJ","SK","SL",
"SM","SN","SO","SP","SQ","SR","SS","ST","SU","SV","SW","SX","SY","SZ","TA","TB",
"TC","TD","TE","TF","TG","TH","TI","TJ","TK","TL","TM","TN","TO","TP","TQ","TR",
"TS","TT","TU","TV","TW","TX","TY","TZ","UA","UB","UC","UD","UE","UF","UG","UH",
"UI","UJ","UK","UL","UM","UN","UO","UP","UQ","UR","US","UT","UU","UV","UW","UX",
"UY","UZ","VA","VB","VC","VD","VE","VF","VG","VH","VI","VJ","VK","VL","VM","VN",
"VO","VP","VQ","VR","VS","VT","VU","VV","VW","VX","VY","VZ","WA","WB","WC","WD",
"WE","WF","WG","WH","WI","WJ","WK","WL","WM","WN","WO","WP","WQ","WR","WS","WT",
"WU","WV","WW","WX","WY","WZ","XA","XB","XC","XD","XE","XF","XG","XH","XI","XJ",
"XK","XL","XM","XN","XO","XP","XQ","XR","XS","XT","XU","XV","XW","XX","XY","XZ",
"YA","YB","YC","YD","YE","YF","YG","YH","YI","YJ","YK","YL","YM","YN","YO","YP",
"YQ","YR","YS","YT","YU","YV","YW","YX","YY","YZ","ZA","ZB","ZC","ZD","ZE","ZF",
"ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV",
"ZW","ZX","ZY","ZZ"
};

/*
    Also tabulate 26 base-26 chars.
*/

/* djb-rwth: removing redundant variables */


/****************************************************************************
 Get a character representing 1st 14-bit triplet
 (bits 0..13 of contiguous array of octets)
****************************************************************************/
const char* base26_triplet_1( const unsigned char *a )
{
    UINT32 b0, b1, h;
    b0 = (UINT32) a[0];             /* 1111 1111  */
#ifndef FIX_BASE26_ENC_BUG
    b1 = (UINT32) ( a[1] & 0x3f );  /* 0011 1111  */
    h = (UINT32) ( b0 | b1 << 8 );
#else
    b1 = (UINT32) ( a[1] & 0xfc );  /* 1111 1100  */
    h = (UINT32) ( b0 << 8 | b1 ) >> 2;
#endif
    return t26[h];
}


/****************************************************************************
 Get a character representing 2nd 14-bit triplet (bits 14..27)
****************************************************************************/
const char* base26_triplet_2( const unsigned char *a )
{
    UINT32 b0, b1, b2, h;
#ifndef FIX_BASE26_ENC_BUG
    b0 = (UINT32) ( a[1] & 0xc0 );   /* 1100 0000  */
    b1 = (UINT32) ( a[2] );         /* 1111 1111  */
    b2 = (UINT32) ( a[3] & 0x0f );  /* 0000 1111  */
    h = (UINT32) ( b0 | b1 << 8 | b2 << 16 ) >> 6;
#else
    b0 = (UINT32) ( a[1] & 0x03 );   /* 0000 0011 */
    b1 = (UINT32) ( a[2] );         /* 1111 1111 */
    b2 = (UINT32) ( a[3] & 0xf0 );  /* 1111 0000 */
    h = (UINT32) ( b0 << 16 | b1 << 8 | b2 ) >> 4;
#endif
    return t26[h];
}


/****************************************************************************
 Get a character representing 3rd 14-bit triplet (bits 28..41)
****************************************************************************/
const char* base26_triplet_3( const unsigned char *a )
{
    UINT32 b0, b1, b2, h;
#ifndef FIX_BASE26_ENC_BUG
    b0 = (UINT32) ( a[3] & 0xf0 );   /* 1111 0000  */
    b1 = (UINT32) ( a[4] );         /* 1111 1111  */
    b2 = (UINT32) ( a[5] & 0x03 );  /* 0000 0011  */
    h = (UINT32) ( b0 | b1 << 8 | b2 << 16 ) >> 4;
#else
    b0 = (UINT32) ( a[3] & 0x0f );   /* 0000 1111 */
    b1 = (UINT32) ( a[4] );         /* 1111 1111  */
    b2 = (UINT32) ( a[5] & 0xc0 );  /* 1100 0000  */
    h = (UINT32) ( b0 << 16 | b1 << 8 | b2 ) >> 6;
#endif
    return t26[h];
}


/****************************************************************************
 Get a character representing 4th 14-bit triplet (bits 42..55)
****************************************************************************/
const char* base26_triplet_4( const unsigned char *a )
{
    UINT32 b0, b1, h;
#ifndef FIX_BASE26_ENC_BUG
    b0 = (UINT32) ( a[5] & 0xfc );   /* 1111 1100  */
    b1 = (UINT32) ( a[6] );         /* 1111 1111  */
    h = (UINT32) ( b0 | b1 << 8 ) >> 2;
#else
    b0 = (UINT32) ( a[5] & 0x3f );   /* 0011 1111  */
    b1 = (UINT32) ( a[6] );         /* 1111 1111  */
    h = (UINT32) ( b0 << 8 | b1 );
#endif
    return t26[h];
}


/*
    Tail dublets
*/

/*
                a4        a3        a2         a1         a0
      28-36:    0001 1111 1111 0000 0000  0000 0000    0000 0000 0000
*/


/****************************************************************************
 Get dublet (bits 28..36)
****************************************************************************/
const char* base26_dublet_for_bits_28_to_36( unsigned char *a )
{
    UINT32 b0, b1, h;
#ifndef FIX_BASE26_ENC_BUG
    b0 = (UINT32) ( a[3] & 0xf0 );    /* 1111 0000  */
    b1 = (UINT32) ( a[4] & 0x1f );    /* 0001 1111  */
    h = (UINT32) ( b0 | b1 << 8 ) >> 4;
#else
    b0 = (UINT32) ( a[3] & 0x0f );    /* 0000 1111  */
    b1 = (UINT32) ( a[4] & 0xf8 );    /* 1111 1000  */
    h = (UINT32) ( b0 << 8 | b1 ) >> 3;
#endif
    return d26[h];
}


/*
                a9        a8        a7
      56-64:    0000 0000 0000 0001 1111 1111
*/


/****************************************************************************
 Get dublet (bits 56..64)
****************************************************************************/
const char* base26_dublet_for_bits_56_to_64( unsigned char *a )
{
    UINT32 b0, b1, h;
#ifndef FIX_BASE26_ENC_BUG
    b0 = (UINT32) ( a[7] );            /* 1111 1111  */
    b1 = (UINT32) ( a[8] & 0x01 );    /* 0000 0001  */
    h = (UINT32) ( b0 | b1 << 8 );
#else
    b0 = (UINT32) ( a[7] );            /* 1111 1111  */
    b1 = (UINT32) ( a[8] & 0x80 );    /* 1000 0000  */
    h = (UINT32) ( b0 << 8 | b1 ) >> 7;
#endif
    return d26[h];
}


/****************************************************************************
 Get hash extension in hexadecimal representation for the major block.
 Len(extension) = 256 - 65 = 191 bit.
****************************************************************************/
void get_xtra_hash_major_hex( const unsigned char *a, char* szXtra )
{
    unsigned char c;
    int i, j, start_byte = 8;
#ifndef FIX_BASE26_ENC_BUG
    c = a[start_byte] & 0xfe;  /* 1111 1110  */
#else
    c = a[start_byte] & 0x7f;  /* 0111 1111  */
#endif
    j = sprintf(szXtra, "%02x", c);
    for (i = start_byte + 1; i < 32; i++)
    {
        j += sprintf(szXtra + j, "%02x", a[i]);
    }
}


/****************************************************************************
 Get hash extension in hexadecimal representation for the minor block.
 Len(extension) = 256 - 37 = 219 bit.
****************************************************************************/
void get_xtra_hash_minor_hex( const unsigned char *a, char* szXtra )
{
    unsigned char c;
    int i, j, start_byte = 4;
#ifndef FIX_BASE26_ENC_BUG
    c = a[start_byte] & 0xe0;  /* 1110 0000  */
#else
    c = a[start_byte] & 0x07;  /* 0000 0111  */
#endif
    j = sprintf(szXtra, "%02x", c);
    for (i = start_byte + 1; i < 32; i++)
    {
        j += sprintf(szXtra + j, "%02x", a[i]);
    }
}
