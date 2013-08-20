#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"
#include "dfibHeap.h"
#include "fibHeap.h"
#include "darray.h"
#include "zlib.h"

#define CNBLOCKSIZE 10000
#define MAXC 10000
#define MAXCinBetween 200

#define MaxNodeInSub 10000
#define GapLowerBound -2000
#define GapUpperBound 300000

#define MaxCntNode 1000

static int DEBUG = 0;
static int DEBUG1 = 0;
static int DEBUG2 = 0;

static int bySmall = 1;

static boolean static_f = 0;
static double OverlapPercent = 0.05;
static double ConflPercent = 0.05;
static int MinWeakCut = 3;
static int gapCounter;
static int orienCounter;
static int orienCounter2;
static int throughCounter;

static int breakPointAtRepeat = 0;
static FILE * snp_fp = NULL;

static DARRAY * solidArray;
static DARRAY * tempArray;

static int solidCounter;

static CTGinHEAP ctg4heapArray[MaxNodeInSub + 1];
static unsigned int nodesInSub[MaxNodeInSub];
static int nodeDistance[MaxNodeInSub];
static int nodeCounter;


static int dCntCounter1;
static int uCntCounter1;
static int dCntCounter2;
static int uCntCounter2;
static unsigned int dCntNodeArr1[MaxCntNode];   //downstream
static unsigned int uCntNodeArr1[MaxCntNode];   //upstream
static int dCntGapArr1[MaxCntNode];
static int uCntGapArr1[MaxCntNode];
static unsigned int dCntNodeArr2[MaxCntNode];
static unsigned int uCntNodeArr2[MaxCntNode];
static int dCntGapArr2[MaxCntNode];
static int uCntGapArr2[MaxCntNode];

static unsigned int * cntNodeArr;
static int * cntGapArr;

static unsigned int nodesInSubInOrder[MaxNodeInSub];
static int nodeDistanceInOrder[MaxNodeInSub];

static DARRAY * scaf3, *scaf5;
static DARRAY * gap3, *gap5;

static unsigned int downstreamCTG[MAXCinBetween];
static unsigned int upstreamCTG[MAXCinBetween];
static int dsCtgCounter;
static int usCtgCounter;

static CONNECT * checkConnect ( unsigned int from_c, unsigned int to_c );
static int maskPuzzle ( int num_connect, unsigned int contigLen );
static void freezing();
static boolean checkOverlapInBetween ( double tolerance );
static int setConnectDelete ( unsigned int from_c, unsigned int to_c, char flag, boolean cleanBinding );
static int setConnectMask ( unsigned int from_c, unsigned int to_c, char mask );
static int setConnectWP ( unsigned int from_c, unsigned int to_c, char flag );

static void general_linearization ( boolean strict );
static void debugging2();
static void smallScaf();
static void clearNewInsFlag();
static void detectBreakScaff();
static void detectBreakScaf();
static boolean checkSimple ( DARRAY * ctgArray, int count );
static void checkCircle();

/*************************************************
Function:
    checkFiles4Scaff
Description:
    Checks the required files for scaffolding.
Input:
    1. infile:      prefix of graph
Output:
    None.
Return:
    1 if all files were OK.
 *************************************************/
boolean checkFiles4Scaff ( char * infile )
{
	char name[7][256];
	boolean filesOK = 1;
	int i = 0;
	sprintf ( name[0], "%s.Arc", infile );
	sprintf ( name[1], "%s.contig", infile );
	sprintf ( name[2], "%s.peGrads", infile );
	sprintf ( name[3], "%s.preGraphBasic", infile );
	sprintf ( name[4], "%s.updated.edge", infile );
	sprintf ( name[5], "%s.readOnContig.gz", infile );

	for ( ; i < 6; i++ )
	{
		filesOK = check_file ( name[i] );

		if ( !filesOK )
		{
			fprintf ( stderr, "%s: no such file or empty file!\n\n", name[i] );
			return filesOK;
		}
	}

	fprintf ( stderr, "Files for scaffold construction are OK.\n\n" );
	return filesOK;
}


/*************************************************
Function:
    getBindCnt
Description:
    Gets the only connection of current contig.
Input:
    1. ctg:     current contig
Output:
    None.
Return:
    The pointer to connection if only one connection was found.
*************************************************/
static CONNECT * getBindCnt ( unsigned int ctg )
{
	CONNECT * ite_cnt;
	CONNECT * bindCnt = NULL;
	CONNECT * temp_cnt = NULL;
	CONNECT * temp3_cnt = NULL;
	int count = 0;
	int count2 = 0;
	int count3 = 0;
	ite_cnt = contig_array[ctg].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->nextInScaf )
		{
			count++;
			bindCnt = ite_cnt;
		}

		if ( ite_cnt->prevInScaf )
		{
			temp_cnt = ite_cnt;
			count2++;
		}

		if ( ite_cnt->singleInScaf )
		{
			temp3_cnt = ite_cnt;
			count3++;
		}

		ite_cnt = ite_cnt->next;
	}

	if ( count == 1 )
		{ return bindCnt; }

	if ( count == 0 && count2 == 1 )
		{ return temp_cnt; }

	if ( count == 0 && count2 == 0 && count3 == 1 )
		{ return temp3_cnt; }

	return NULL;
}

/*************************************************
Function:
    createAnalogousCnt
Description:
    Creats new relation between two non-connected contigs
    according to toplogy sturcture supported by large insert size
    paired-end reads.
Input:
    1. sourceStart: contig that has connection relations to other two
                        contigs
    2. originCnt:       direct connection supported by large insert size
                            paired-end reads
    3. gap:         calculated distance between two non-connected contigs
    4. targetStart:     left contig of those two non-connected contigs
    5. targetStop:      right contig of those two non-connected contigs
Output:
    None.
Return:
    1 if creation successed.
*************************************************/
static boolean createAnalogousCnt ( unsigned int sourceStart,
                                    CONNECT * originCnt, int gap,
                                    unsigned int targetStart, unsigned int targetStop )
{
	CONNECT * temp_cnt;
	unsigned int balTargetStart = getTwinCtg ( targetStart );
	unsigned int balTargetStop = getTwinCtg ( targetStop );
	unsigned int balSourceStart = getTwinCtg ( sourceStart );
	unsigned int balSourceStop = getTwinCtg ( originCnt->contigID );
	boolean change_flag = 0;
	int add_weight = originCnt->weight;

	if ( gap < GapLowerBound )
	{
		gapCounter++;
		originCnt->deleted = 1;
		temp_cnt = getCntBetween ( balSourceStop, balSourceStart );
		temp_cnt->deleted = 1;
		return change_flag;
	}

	int startLen = ( int ) contig_array[targetStart].length;
	int stopLen = ( int ) contig_array[targetStop].length;

	if ( gap < -overlaplen )
	{
		if ( ( int ) contig_array[targetStart].length - overlaplen <= -gap && ( int ) contig_array[targetStart].length > 5 * overlaplen )
		{
			unsigned int tmp_id = targetStart;
			targetStart = targetStop;
			targetStop = tmp_id;
			tmp_id = balTargetStart;
			balTargetStart = balTargetStop;
			balTargetStop = tmp_id;
			gap = -gap;
			change_flag = 1;
		}
		else
		{
			gapCounter++;
			return change_flag;
		}
	}

	originCnt->deleted = 1;
	temp_cnt = getCntBetween ( balSourceStop, balSourceStart );
	temp_cnt->deleted = 1;
	temp_cnt = add1Connect ( targetStart, targetStop, gap, add_weight, 1 );

	if ( temp_cnt )
		{ temp_cnt->inherit = 1; }

	temp_cnt = add1Connect ( balTargetStop, balTargetStart, gap, add_weight, 1 );

	if ( temp_cnt )
		{ temp_cnt->inherit = 1; }

	return change_flag;
}


/*************************************************
Function:
    add1LongPEcov
Description:
    Increases the connections weight supported by large insert size
    paired-end reads between contigs in one scaffold.
Input:
    1. fromCtg:     left contig of pair contig connected by large insert
                            size paired-end reads
    2. toCtg:           right contig of pair contig connected by large insert
                            size paired-end reads
    3. weight:      weight of connection
Output:
    None.
Return:
    None.
*************************************************/
static void add1LongPEcov ( unsigned int fromCtg, unsigned int toCtg, int weight )
{
	//check if they are on the same scaff
	if ( contig_array[fromCtg].from_vt != contig_array[toCtg].from_vt ||
	        contig_array[fromCtg].to_vt != contig_array[toCtg].to_vt )
	{
		fprintf ( stderr, "Warning from add1LongPEcov: contig %d and %d not on the same scaffold\n",
		          fromCtg, toCtg );
		return;
	}

	if ( contig_array[fromCtg].indexInScaf >= contig_array[toCtg].indexInScaf )
	{
		fprintf ( stderr, "Warning from add1LongPEcov: wrong about order between contig %d and %d\n",
		          fromCtg, toCtg );
		return;
	}

	CONNECT * bindCnt;
	unsigned int prevCtg = fromCtg;
	bindCnt = getBindCnt ( fromCtg );

	while ( bindCnt )
	{
		if ( bindCnt->maxGap + weight <= 1000 )
			{ bindCnt->maxGap += weight; }
		else
			{ bindCnt->maxGap = 1000; }

		if ( fromCtg == 0 && toCtg == 0 )
			fprintf ( stderr, "link (%d %d ) covered by link (%d %d), wt %d\n",
			          prevCtg, bindCnt->contigID, fromCtg, toCtg, weight );

		if ( bindCnt->contigID == toCtg )
			{ break; }

		prevCtg = bindCnt->contigID;
		bindCnt = bindCnt->nextInScaf;
	}

	unsigned int bal_fc = getTwinCtg ( fromCtg );
	unsigned int bal_tc = getTwinCtg ( toCtg );
	bindCnt = getBindCnt ( bal_tc );
	prevCtg = bal_tc;

	while ( bindCnt )
	{
		if ( bindCnt->maxGap + weight <= 1000 )
			{ bindCnt->maxGap += weight; }
		else
			{ bindCnt->maxGap = 1000; }

		if ( fromCtg == 0 && toCtg == 0 )
			fprintf ( stderr, "link (%d %d ) covered by link (%d %d), wt %d\n",
			          prevCtg, bindCnt->contigID, fromCtg, toCtg, weight );

		if ( bindCnt->contigID == bal_fc )
			{ return; }

		prevCtg = bindCnt->contigID;
		bindCnt = bindCnt->nextInScaf;
	}
}

/*************************************************
Function:
    downSlide
Description:
    Deals with connections supported by large insert size paired-end
    reads. If a connecton is abnormal, or the two connected contigs
    are in the same scaffold and the distance between them is normal,
    set the connection's status to 'deleted'. If the two connected
    contigs locate in different scaffolds, trys to creat a connection
    between these two scaffolds.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void downSlide()
{
	int len = 0, gap;
	unsigned int i;
	CONNECT * ite_cnt, *bindCnt, *temp_cnt;
	unsigned int bottomCtg, topCtg, bal_i;
	unsigned int targetCtg, bal_target;
	boolean getThrough, orienConflict;
	int slideLen, slideLen2;
	int slideexchange1 = 0, slideexchange2 = 0;
	orienCounter = throughCounter = 0;
	orienCounter2 = 0;
	int slidebreak1 = 0, slidebreak2 = 0, slidebreak = 0, recoverCnt = 0;

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].mask || !contig_array[i].downwardConnect )
			{ continue; }

		bindCnt = getBindCnt ( i );

		if ( !bindCnt )
			{ continue; }

		bal_i = getTwinCtg ( i );
		len = slideLen = 0;
		bottomCtg = i;

		//find the last unmasked contig in this binding
		while ( bindCnt->nextInScaf )
		{
			len += bindCnt->gapLen + contig_array[bindCnt->contigID].length;

			if ( contig_array[bindCnt->contigID].mask == 0 )
			{
				bottomCtg = bindCnt->contigID;
				slideLen = len;
			}

			bindCnt = bindCnt->nextInScaf;
		}

		len += bindCnt->gapLen + contig_array[bindCnt->contigID].length;

		if ( contig_array[bindCnt->contigID].mask == 0 || bottomCtg == 0 )
		{
			bottomCtg = bindCnt->contigID;
			slideLen = len;
		}

		//check each connetion from long pair ends
		ite_cnt = contig_array[i].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->deleted || ite_cnt->mask || ite_cnt->singleInScaf
			        || ite_cnt->nextInScaf || ite_cnt->prevInScaf || ite_cnt->inherit )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			targetCtg = ite_cnt->contigID;

			if ( contig_array[i].from_vt == contig_array[targetCtg].from_vt ) // on the same scaff
			{
				if ( contig_array[i].indexInScaf > contig_array[targetCtg].indexInScaf )
					{ orienCounter++; }
				else
					{ throughCounter++; }

				setConnectDelete ( i, ite_cnt->contigID, 1, 0 );
				ite_cnt = ite_cnt->next;
				continue;
			}

			if ( ( ins_var_idx > 0 ) && ( slideLen > Insert_size * ins_var_idx ) )
			{
				setConnectDelete ( i, ite_cnt->contigID, 1, 0 );
				ite_cnt = ite_cnt->next;
				slidebreak1++;
				continue;
			}

			//contig i and targetctg is not in same scaffold
			//check if this connection conflicts with previous scaffold orientationally
			temp_cnt = getBindCnt ( targetCtg );
			orienConflict = 0;

			if ( temp_cnt )
			{
				while ( temp_cnt->nextInScaf )
				{
					if ( temp_cnt->contigID == i )
					{
						orienConflict = 1;
						fprintf ( stderr, "Warning from downSlide: still on the same scaff: %d and %d\n"
						          , i, targetCtg );
						fprintf ( stderr, "on scaff %d and %d\n",
						          contig_array[i].from_vt, contig_array[targetCtg].from_vt );
						fprintf ( stderr, "on bal_scaff %d and %d\n",
						          contig_array[bal_target].to_vt, contig_array[bal_i].to_vt );
						break;
					}

					temp_cnt = temp_cnt->nextInScaf;
				}

				if ( temp_cnt->contigID == i )
					{ orienConflict = 1; }
			}

			if ( orienConflict )
			{
				orienCounter++;
				orienCounter2++;
				setConnectDelete ( i, ite_cnt->contigID, 1, 0 );
				ite_cnt = ite_cnt->next;
				continue;
			}

			//connection path to i was not found
			//find the most top contig along previous scaffold starting with the target contig of this connection
			bal_target = getTwinCtg ( targetCtg );
			slideLen2 = 0;

			if ( contig_array[targetCtg].mask == 0 )
			{
				topCtg = bal_target;
			}
			else
			{
				topCtg = 0;
			}

			temp_cnt = getBindCnt ( bal_target );
			getThrough = len = 0;
			int slidebreak = 0;

			if ( temp_cnt )
			{
				//find the last contig in this binding
				while ( temp_cnt->nextInScaf )
				{
					//check if this route reaches bal_i
					if ( temp_cnt->contigID == bal_i )
					{
						fprintf ( stderr, "Warning from downSlide: (B) still on the same scaff: %d and %d (%d and %d)\n",
						          i, targetCtg, bal_target, bal_i );
						fprintf ( stderr, "on scaff %d and %d\n",
						          contig_array[i].from_vt, contig_array[targetCtg].from_vt );
						fprintf ( stderr, "on bal_scaff %d and %d\n",
						          contig_array[bal_target].to_vt, contig_array[bal_i].to_vt );
						getThrough = 1;
						break;
					}

					len += temp_cnt->gapLen + contig_array[temp_cnt->contigID].length;

					if ( contig_array[temp_cnt->contigID].mask == 0 )
					{
						topCtg = temp_cnt->contigID;
						slideLen2 = len;
					}

					if ( ( ins_var_idx > 0 ) && ( len > ins_var_idx * Insert_size ) )
					{
						slidebreak = 1;
						break;
					}

					temp_cnt = temp_cnt->nextInScaf;
				}

				len += temp_cnt->gapLen + contig_array[temp_cnt->contigID].length;

				if ( contig_array[temp_cnt->contigID].mask == 0 || topCtg == 0 )
				{
					topCtg = temp_cnt->contigID;
					slideLen2 = len;
				}

				if ( slidebreak == 1 )
				{
					setConnectDelete ( i, ite_cnt->contigID, 1, 0 );
					ite_cnt = ite_cnt->next;
					slidebreak2++;
					continue;
				}

				if ( temp_cnt->contigID == bal_i )
					{ getThrough = 1; }
				else
					{ topCtg = getTwinCtg ( topCtg ); }
			}
			else
				{ topCtg = targetCtg; }

			if ( getThrough )
			{
				throughCounter++;
				setConnectDelete ( i, ite_cnt->contigID, 1, 0 );
				ite_cnt = ite_cnt->next;
				continue;
			}

			//connection path to bal_id was not found
			CONNECT * dh_cnt;
			gap = ite_cnt->gapLen - slideLen - slideLen2;
			dh_cnt = getCntBetween ( topCtg, bottomCtg );

			if ( dh_cnt && dh_cnt->weight >= MinWeakCut )
			{
				slideexchange1++;
				setConnectDelete ( topCtg, bottomCtg, 0, 0 );
				setConnectMask ( topCtg, bottomCtg, 0 );
				ite_cnt = ite_cnt->next;
				continue;
			}

			//add a connection between bottomCtg and topCtg
			if ( bottomCtg != topCtg && ! ( i == bottomCtg && targetCtg == topCtg ) )
			{
				boolean creat_flag = createAnalogousCnt ( i, ite_cnt, gap, bottomCtg, topCtg );

				if ( creat_flag )
					{ slideexchange2++; }

				if ( contig_array[bottomCtg].mask || contig_array[topCtg].mask )
					{ fprintf ( stderr, "downSlide to masked contig, bottomCtg %u[mask %d], topCtg %u[mask %d]\n", bottomCtg, contig_array[bottomCtg].mask, topCtg, contig_array[topCtg].mask ); }
			}

			ite_cnt = ite_cnt->next;
		}
	}

	//    fprintf(stderr,"downSliding stat:\norienConflict\tfall_inside\tslidebreak1\tslidebreak2\trecoverCnt\tslideexchange1\tslideexchange2\n%d\t%d\t%d\t%d\t%d\t%d\t%d\n",orienCounter, throughCounter, slidebreak1, slidebreak2, recoverCnt,  slideexchange1, slideexchange2);
	fprintf ( stderr, "Add large insert size PE links: %d orientation-conflict links, %d contigs acrossed by normal links.\n", orienCounter, throughCounter );
}

/*************************************************
Function:
    setNextInScaf
Description:
    Sets the downstream connection of current connection.
Input:
    1. cnt:         current connection
    2. nextCnt:     next connection
Output:
    None.
Return:
    1 if setting successed.
*************************************************/
static boolean setNextInScaf ( CONNECT * cnt, CONNECT * nextCnt )
{
	if ( !cnt )
	{
		fprintf ( stderr, "setNextInScaf: empty pointer\n" );
		return 0;
	}

	if ( !nextCnt )
	{
		cnt->nextInScaf = nextCnt;
		return 1;
	}

	if ( cnt->mask || cnt->deleted )
	{
		fprintf ( stderr, "setNextInScaf: cnt is masked or deleted\n" );
		return 0;
	}

	if ( nextCnt->deleted || nextCnt->mask )
	{
		fprintf ( stderr, "setNextInScaf: nextCnt is masked or deleted\n" );
		return 0;
	}

	cnt->nextInScaf = nextCnt;
	return 1;
}

/*************************************************
Function:
    setPrevInScaf
Description:
    Sets the upstream connection status of current connection.
Input:
    1. cnt:     current connection
    2. flag:        new status
Output:
    None.
Return:
    1 if setting successed.
*************************************************/
static boolean setPrevInScaf ( CONNECT * cnt, boolean flag )
{
	if ( !cnt )
	{
		fprintf ( stderr, "setPrevInScaf: empty pointer\n" );
		return 0;
	}

	if ( !flag )
	{
		cnt->prevInScaf = flag;
		return 1;
	}

	if ( cnt->mask || cnt->deleted )
	{
		fprintf ( stderr, "setPrevInScaf: cnt is masked or deleted\n" );
		return 0;
	}

	cnt->prevInScaf = flag;
	return 1;
}


/*************************************************
Function:
    substitueUSinScaf
Description:
    Substitutes the upstream connection of current connection with
    new connection.
        from_c
                        -> branch_c -> to_c
            from_c_new
Input:
    1. origin:          original upstream connection of current connection
                            (from_c -> branch_c)
    2. from_c_new:  new upstream contig of current contig (branch_c)
Output:
    None.
Return:
    None.
*************************************************/
static void substitueUSinScaf ( CONNECT * origin, unsigned int from_c_new )
{
	if ( !origin || !origin->nextInScaf )
		{ return; }

	unsigned int branch_c, to_c;
	unsigned int bal_branch_c, bal_to_c;
	unsigned int bal_from_c_new = getTwinCtg ( from_c_new );
	CONNECT * bal_origin, *bal_nextCNT, *prevCNT, *bal_prevCNT;
	branch_c = origin->contigID;
	to_c = origin->nextInScaf->contigID;
	bal_branch_c = getTwinCtg ( branch_c );
	bal_to_c = getTwinCtg ( to_c );
	prevCNT = checkConnect ( from_c_new, branch_c );
	bal_nextCNT = checkConnect ( bal_to_c, bal_branch_c );

	if ( !bal_nextCNT )
	{
		fprintf ( stderr, "substitueUSinScaf: no connect between %d and %d\n", bal_to_c, bal_branch_c );
		return;
	}

	bal_origin = bal_nextCNT->nextInScaf;
	bal_prevCNT = checkConnect ( bal_branch_c, bal_from_c_new );
	setPrevInScaf ( bal_nextCNT->nextInScaf, 0 );
	setNextInScaf ( prevCNT, origin->nextInScaf );
	setNextInScaf ( bal_nextCNT, bal_prevCNT );
	setPrevInScaf ( bal_prevCNT, 1 );
	setNextInScaf ( origin, NULL );
	setPrevInScaf ( bal_origin, 0 );
}


/*************************************************
Function:
    substitueDSinScaf
Description:
    Substitutes the downstream connection of current connection with
    new connection.
                                    to_c
            from_c -> branch_c ->
                                  to_c_new
Input:
    1. origin:          original downstream connection of current connection
                         (branch_c -> to_c)
    2. branch_c:        current contig
    3. to_c_new:        new downstream contig of current contig
Output:
    None.
Return:
    None.
*************************************************/
static void substitueDSinScaf ( CONNECT * origin, unsigned int branch_c, unsigned int to_c_new )
{
	if ( !origin || !origin->prevInScaf )
		{ return; }

	unsigned int to_c;
	unsigned int bal_branch_c, bal_to_c, bal_to_c_new;
	unsigned int from_c, bal_from_c;
	CONNECT * bal_origin, *prevCNT, *bal_prevCNT;
	CONNECT * nextCNT, *bal_nextCNT;
	to_c = origin->contigID;
	bal_branch_c = getTwinCtg ( branch_c );
	bal_to_c = getTwinCtg ( to_c );
	bal_origin = getCntBetween ( bal_to_c, bal_branch_c );

	if ( !bal_origin )
	{
		fprintf ( stderr, "substitueDSinScaf: no connect between %d and %d\n", bal_to_c, bal_branch_c );
		return;
	}

	if ( bal_origin->nextInScaf )
		{ bal_from_c = bal_origin->nextInScaf->contigID; }
	else
	{
		fprintf ( stderr, "next null! %d\t%d\n", bal_to_c, bal_branch_c );
		exit ( 3 );
	}

	bal_from_c = bal_origin->nextInScaf->contigID;
	from_c = getTwinCtg ( bal_from_c );
	bal_to_c_new = getTwinCtg ( to_c_new );
	prevCNT = checkConnect ( from_c, branch_c );
	nextCNT = checkConnect ( branch_c, to_c_new );
	setNextInScaf ( prevCNT, nextCNT );
	setPrevInScaf ( nextCNT, 1 );
	bal_nextCNT = checkConnect ( bal_to_c_new, bal_branch_c );
	bal_prevCNT = checkConnect ( bal_branch_c, bal_from_c );
	setNextInScaf ( bal_nextCNT, bal_prevCNT );
	setPrevInScaf ( origin, 0 );
	setNextInScaf ( bal_origin, NULL );
}

/*************************************************
Function:
    validConnect
Description:
    Calculates the valid connections number of a contig.
    1. Check whether the contig has upstream conntcion and
        downstream connection.
    2. Calculates the non-deleted and non-masked connections number.
Input:
    1. ctg:         contig
    2. preCNT:      upstream connection
Output:
    None.
Return:
    0 if contig did NOT have connection to other contigs.
    1 if contig had upstream conntcion and downstream connection.
    Non-deleted and non-masked connections number otherwise.
*************************************************/
static int validConnect ( unsigned int ctg, CONNECT * preCNT )
{
	if ( preCNT && preCNT->nextInScaf )
		{ return 1; }

	CONNECT * cn_temp;
	int count = 0;

	if ( !contig_array[ctg].downwardConnect )
		{ return count; }

	cn_temp = contig_array[ctg].downwardConnect;

	while ( cn_temp )
	{
		if ( !cn_temp->deleted && !cn_temp->mask )
			{ count++; }

		cn_temp = cn_temp->next;
	}

	return count;
}

/*************************************************
Function:
    getNextContig
Description:
    Gets the connection to next contig of current contig through constructed
    connection in   scaffold or direct connection to other contigs.
    In the latter case, one and only one non-deleted and non-masked
    connection is qualified.
Input:
    1. ctg:         current contig
    2. preCNT:      upstream connection of current contig
Output:
    1. exception:       indicate that whether found connection through
        constructed connection in scaffold is deleted or masked
Return:
    Pointer to qualified connection or NULL.
*************************************************/
static CONNECT * getNextContig ( unsigned int ctg, CONNECT * preCNT, boolean * exception )
{
	CONNECT * cn_temp, *retCNT = NULL, *dh_cnt;
	int count = 0, valid_in;
	unsigned int nextCtg, bal_ctg;
	*exception = 0;

	if ( preCNT && preCNT->nextInScaf )
	{
		if ( preCNT->contigID != ctg )
			{ fprintf ( stderr, "pre cnt does not lead to %d\n", ctg ); }

		nextCtg = preCNT->nextInScaf->contigID;
		cn_temp = getCntBetween ( ctg, nextCtg );
		dh_cnt = getCntBetween ( getTwinCtg ( nextCtg ), getTwinCtg ( ctg ) );

		if ( cn_temp && ( cn_temp->mask || cn_temp->deleted ) )
		{
			int id1 = 0, id2 = 0;

			if ( dh_cnt->nextInScaf )
			{
				id1 = dh_cnt->nextInScaf->contigID;
				id2 = getTwinCtg ( dh_cnt->nextInScaf->contigID );
			}

			if ( !cn_temp->prevInScaf )
				{ fprintf ( stderr, "not even has a prevInScaf %d and %d, %d and %d with before %d with twin %d\n", ctg, nextCtg, getTwinCtg ( nextCtg ), getTwinCtg ( ctg ), id1, id2 ); }

			cn_temp = getCntBetween ( getTwinCtg ( nextCtg ),
			                          getTwinCtg ( ctg ) );

			if ( !cn_temp->nextInScaf )
				{ fprintf ( stderr, "its twin cnt not  has a nextInScaf\n" ); }

			fflush ( stdout );
			*exception = 1;
		}
		else
			{ return preCNT->nextInScaf; }
	}

	bal_ctg = getTwinCtg ( ctg );
	valid_in = validConnect ( bal_ctg, NULL );

	if ( valid_in > 1 )
		{ return NULL; }

	if ( !contig_array[ctg].downwardConnect )
		{ return NULL; }

	cn_temp = contig_array[ctg].downwardConnect;

	while ( cn_temp )
	{
		if ( cn_temp->mask || cn_temp->deleted )
		{
			cn_temp = cn_temp->next;
			continue;
		}

		count++;

		if ( count == 1 )
			{ retCNT = cn_temp; }
		else if ( count == 2 )
			{ return NULL; }

		cn_temp = cn_temp->next;
	}

	return retCNT;
}

/*************************************************
Function:
    checkConnect
Description:
    Check whether the connection between two contigs exists and
    the connection's status.
Input:
    1. from_c:      left contig
    2. to_c:            right contig
Output:
    None.
Return:
    Pointer to qualified connection or NULL.
*************************************************/
static CONNECT * checkConnect ( unsigned int from_c, unsigned int to_c )
{
	CONNECT * cn_temp = getCntBetween ( from_c, to_c );

	if ( !cn_temp )
		{ return NULL; }

	if ( !cn_temp->mask && !cn_temp->deleted )
		{ return cn_temp; }

	//else
	//printf("masked or deleted: %d\t%d\t%d\t%d\n",from_c,to_c,cn_temp->mask,cn_temp->deleted);
	return NULL;
}

/*************************************************
Function:
    setConnectMask
Description:
    Sets the "mask" status of connection between two contigs and
    releated upstream and downstream connection's status if these
    connections exist.
Input:
    1. from_c:      left contig
    2. to_c:            right contig
    3. mask:            new status
Output:
    None.
Return:
    1 if setting successed.
*************************************************/
static int setConnectMask ( unsigned int from_c, unsigned int to_c, char mask )
{
	CONNECT * cn_temp, *cn_bal, *cn_ds, *cn_us;
	unsigned int bal_fc = getTwinCtg ( from_c );
	unsigned int bal_tc = getTwinCtg ( to_c );
	unsigned int ctg3, bal_ctg3;
	cn_temp = getCntBetween ( from_c, to_c );
	cn_bal = getCntBetween ( bal_tc, bal_fc );

	if ( !cn_temp || !cn_bal )
	{
		return 0;
	}

	cn_temp->mask = mask;
	cn_bal->mask = mask;

	if ( !mask )
		{ return 1; }

	if ( cn_temp->nextInScaf ) //undo the binding
	{
		setPrevInScaf ( cn_temp->nextInScaf, 0 );
		ctg3 = cn_temp->nextInScaf->contigID;
		setNextInScaf ( cn_temp, NULL );
		bal_ctg3 = getTwinCtg ( ctg3 );
		cn_ds = getCntBetween ( bal_ctg3, bal_tc );
		setNextInScaf ( cn_ds, NULL );
		setPrevInScaf ( cn_bal, 0 );
	}

	// ctg3 -> from_c -> to_c
	// bal_ctg3 <- bal_fc <- bal_tc
	if ( cn_bal->nextInScaf )
	{
		setPrevInScaf ( cn_bal->nextInScaf, 0 );
		bal_ctg3 = cn_bal->nextInScaf->contigID;
		setNextInScaf ( cn_bal, NULL );
		ctg3 = getTwinCtg ( bal_ctg3 );
		cn_us = getCntBetween ( ctg3, from_c );
		setNextInScaf ( cn_us, NULL );
		setPrevInScaf ( cn_temp, 0 );
	}

	return 1;
}

/*************************************************
Function:
    setConnectUsed
Description:
    Sets the "used" status of connection between two contigs if the
    connection exists.
Input:
    1. from_c:      left contig
    2. to_c:            right contig
    3. flag:            new status
Output:
    None.
Return:
    1 if setting successed.
*************************************************/
static boolean setConnectUsed ( unsigned int from_c, unsigned int to_c, char flag )
{
	CONNECT * cn_temp, *cn_bal;
	unsigned int bal_fc = getTwinCtg ( from_c );
	unsigned int bal_tc = getTwinCtg ( to_c );
	cn_temp = getCntBetween ( from_c, to_c );
	cn_bal = getCntBetween ( bal_tc, bal_fc );

	if ( !cn_temp || !cn_bal )
	{
		return 0;
	}

	cn_temp->used = flag;
	cn_bal->used = flag;
	return 1;
}

/*************************************************
Function:
    setConnectWP
Description:
    Sets the "weakPoint" status of connection between two contigs if the
    connection exists.
Input:
    1. from_c:      left contig
    2. to_c:            right contig
    3. flag:            new status
Output:
    None.
Return:
    1 if setting successed.
*************************************************/
static int setConnectWP ( unsigned int from_c, unsigned int to_c, char flag )
{
	CONNECT * cn_temp, *cn_bal;
	unsigned int bal_fc = getTwinCtg ( from_c );
	unsigned int bal_tc = getTwinCtg ( to_c );
	cn_temp = getCntBetween ( from_c, to_c );
	cn_bal = getCntBetween ( bal_tc, bal_fc );

	if ( !cn_temp || !cn_bal )
	{
		return 0;
	}

	cn_temp->weakPoint = flag;
	cn_bal->weakPoint = flag;
	return 1;
}

/*************************************************
Function:
    setConnectDelete
Description:
    Sets the "deleted" status of connection between two contigs if the
    connection exists. Cleans the related upstream and downstream
    connections if required.
Input:
    1. from_c:          left contig
    2. to_c:                right contig
    3. flag:                new status
    4. cleanBinding:        indicator of whether cleaning related connections
Output:
    None.
Return:
    1 if setting successed.
*************************************************/
static int setConnectDelete ( unsigned int from_c, unsigned int to_c, char flag, boolean cleanBinding )
{
	CONNECT * cn_temp, *cn_bal;
	unsigned int bal_fc = getTwinCtg ( from_c );
	unsigned int bal_tc = getTwinCtg ( to_c );
	cn_temp = getCntBetween ( from_c, to_c );
	cn_bal = getCntBetween ( bal_tc, bal_fc );

	if ( !cn_temp || !cn_bal )
	{
		return 0;
	}

	cn_temp->deleted = flag;
	cn_bal->deleted = flag;

	if ( !flag )
		{ return 1; }

	if ( cleanBinding )
	{
		cn_temp->prevInScaf = 0;
		cn_temp->nextInScaf = NULL;
		cn_bal->prevInScaf = 0;
		cn_bal->nextInScaf = NULL;
	}

	return 1;
}

/*************************************************
Function:
    maskContig
Description:
    Sets the "mask" status of target contig and its connections to
    other contigs except those connections in scaffold.
Input:
    1. ctg:     target contig
    2. flag:        new status
Output:
    None.
Return:
    None.
*************************************************/
static void maskContig ( unsigned int ctg, boolean flag )
{
	unsigned int bal_ctg, ctg2, bal_ctg2;
	CONNECT * cn_temp;
	bal_ctg = getTwinCtg ( ctg );
	cn_temp = contig_array[ctg].downwardConnect;

	while ( cn_temp )
	{
		if ( cn_temp->mask || cn_temp->prevInScaf || cn_temp->nextInScaf || cn_temp->singleInScaf )
		{
			cn_temp = cn_temp->next;
			continue;
		}

		ctg2 = cn_temp->contigID;
		setConnectMask ( ctg, ctg2, flag );
		cn_temp = cn_temp->next;
	}

	// bal_ctg2 <- bal_ctg
	cn_temp = contig_array[bal_ctg].downwardConnect;

	while ( cn_temp )
	{
		if ( cn_temp->mask || cn_temp->prevInScaf || cn_temp->nextInScaf || cn_temp->singleInScaf )
		{
			cn_temp = cn_temp->next;
			continue;
		}

		bal_ctg2 = cn_temp->contigID;
		setConnectMask ( bal_ctg, bal_ctg2, flag );
		cn_temp = cn_temp->next;
	}

	contig_array[ctg].mask = flag;
	contig_array[bal_ctg].mask = flag;
}


/*************************************************
Function:
    maskPuzzle
Description:
    For contigs longer than "contigLen", masks them if they have
    more than one valid connection in either upstream of downstream
    direction, and have more than 'num_connect' valid connection.
Input:
    1. num_connect:     allowed valid connection number cutoff
    2. contigLen:           contig length cutoff
Output:
    None.
Return:
    Masked contig number.
*************************************************/
static int maskPuzzle ( int num_connect, unsigned int contigLen )
{
	int in_num, out_num, flag = 0, puzzleCounter = 0;
	unsigned int i, bal_i;
	fprintf ( stderr, "Start to mask puzzles.\n" );

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contigLen && contig_array[i].length > contigLen )
			{ break; }

		if ( contig_array[i].mask )
			{ continue; }

		bal_i = getTwinCtg ( i );
		in_num = validConnect ( bal_i, NULL );
		out_num = validConnect ( i, NULL );

		if ( ( in_num > 1 || out_num > 1 ) && ( in_num + out_num >= num_connect ) )
		{
			flag++;
			maskContig ( i, 1 );
		}

		// upstream connection in scaffold
		in_num = validConnect ( bal_i, NULL );
		// downstream connection in scaffold
		out_num = validConnect ( i, NULL );

		if ( in_num > 1 || out_num > 1 )
		{
			puzzleCounter++;
			//debugging2(i);
		}

		if ( isSmallerThanTwin ( i ) )
			{ i++; }
	}

	fprintf ( stderr, " Masked contigs      %d\n Remained puzzles    %d\n", flag, puzzleCounter );
	return flag;
}

/*************************************************
Function:
    deleteWeakCnt
Description:
    Updates the status of contigs' connections according to "cut_off",
    and then mask those contigs which form circle structure.

    A -> B -> A

Input:
    1. cut_off:     weight cutoff of connection
Output:
    None.
Return:
    None.
*************************************************/
static void deleteWeakCnt ( int cut_off )
{
	unsigned int i;
	CONNECT * cn_temp1;
	int weaks = 0, counter = 0;

	for ( i = 1; i <= num_ctg; i++ )
	{
		cn_temp1 = contig_array[i].downwardConnect;

		while ( cn_temp1 )
		{
			if ( !cn_temp1->mask && !cn_temp1->deleted && !cn_temp1->nextInScaf
			        && !cn_temp1->singleInScaf && !cn_temp1->prevInScaf )
			{
				counter++;
			}

			if ( cn_temp1->weak && cn_temp1->deleted && cn_temp1->weight >= cut_off )
			{
				cn_temp1->deleted = 0;
				cn_temp1->weak = 0;
			}
			else if ( !cn_temp1->deleted && cn_temp1->weight > 0 && cn_temp1->weight < cut_off
			          && !cn_temp1->nextInScaf && !cn_temp1->prevInScaf )
			{
				cn_temp1->deleted = 1;
				cn_temp1->weak = 1;

				if ( cn_temp1->singleInScaf )
					{ cn_temp1->singleInScaf = 0; }

				if ( !cn_temp1->mask )
					{ weaks++; }
			}

			cn_temp1 = cn_temp1->next;
		}
	}

	if ( counter > 0 )
		{ fprintf ( stderr, " Active connections    %d\n Weak connections      %d\n Weak ratio            %.1f%%\n", counter, weaks, ( float ) weaks / counter * 100 ); }

	checkCircle();
}


/*************************************************
Function:
    linearC2C
Description:
    For three contigs A, B and C, connections A->B and A->C exist,
    check if there is a connection path connecting B and C or
    whether a new connection can be created between B and C.

       ------->
    A -> B -?-> C

Input:
    1. starter:     leftmost contig
    2. cnt2c1:      connection A->B
    3. c2:          rightmost contig
    4. min_dis:     minimum distance between middle and rightmost contig
    5. max_dis:     maximum distance between middle and rightmost contig
Output:
    None.
Return:
    -1 if B and C were the same contig.
    1 if connection path was found or connection was created.
    0 otherwise.
*************************************************/
static int linearC2C ( unsigned int starter, CONNECT * cnt2c1, unsigned int c2, int min_dis, int max_dis )
{
	int out_num, in_num;
	CONNECT * prevCNT, *cnt, *cn_temp;
	unsigned int c1, bal_c1, ctg, bal_c2;
	int len = 0;
	unsigned int bal_start = getTwinCtg ( starter );
	boolean excep;
	c1 = cnt2c1->contigID;

	if ( c1 == c2 )
	{
		fprintf ( stderr, "linearC2C: c1(%d) and c2(%d) are the same contig.\n", c1, c2 );
		return -1;
	}

	bal_c1 = getTwinCtg ( c1 );
	dsCtgCounter = 1;
	usCtgCounter = 0;
	downstreamCTG[dsCtgCounter++] = c1;
	bal_c2 = getTwinCtg ( c2 );
	upstreamCTG[usCtgCounter++] = bal_c2;
	// check if c1 is linearly connected to c2 by pe connections
	cnt = prevCNT = cnt2c1;

	while ( ( cnt = getNextContig ( c1, prevCNT, &excep ) ) != NULL )
	{
		c1 = cnt->contigID;
		len += cnt->gapLen + contig_array[c1].length;

		if ( c1 == c2 )
		{
			usCtgCounter--;
			return 1; //is interleaving.
		}

		if ( len > max_dis || c1 == starter || c1 == bal_start )
			{ return 0; }

		downstreamCTG[dsCtgCounter++] = c1;

		if ( dsCtgCounter >= MAXCinBetween )
		{
			fprintf ( stderr, "%d downstream contigs, start at %d, max_dis %d, current dis %d\n"
			          , dsCtgCounter, starter, max_dis, len );
			return 0;
		}

		prevCNT = cnt;
	}

	out_num = validConnect ( c1, NULL ); //new c1 should have no outgoing link.

	if ( out_num )
		{ return 0; }

	//find the most upstream contig to c2
	cnt = prevCNT = NULL;
	ctg = bal_c2;

	while ( ( cnt = getNextContig ( ctg, prevCNT, &excep ) ) != NULL )
	{
		ctg = cnt->contigID;
		len += cnt->gapLen + contig_array[ctg].length;

		if ( len > max_dis || ctg == starter || ctg == bal_start )
			{ return 0; }

		prevCNT = cnt;
		upstreamCTG[usCtgCounter++] = ctg;

		if ( usCtgCounter >= MAXCinBetween )
		{
			fprintf ( stderr, "%d upstream contigs, start at %d, max_dis %d, current dis %d\n"
			          , usCtgCounter, starter, max_dis, len );
			return 0;
		}
	}

	if ( dsCtgCounter + usCtgCounter > MAXCinBetween )
	{
		fprintf ( stderr, "%d downstream and %d upstream contigs.\n", dsCtgCounter, usCtgCounter );
		return 0;
	}

	out_num = validConnect ( ctg, NULL ); //new c2 have no incoming link.

	if ( out_num )
	{
		return 0;
	}

	c2 = getTwinCtg ( ctg );
	min_dis -= len;
	max_dis -= len;

	if ( c1 == c2 || c1 == ctg || max_dis < 0 )
		{ return 0; }

	usCtgCounter--;
	cn_temp = getCntBetween ( c1, c2 ); //have connection between new c1 and new c2

	if ( cn_temp )
	{
		setConnectMask ( c1, c2, 0 );
		setConnectDelete ( c1, c2, 0, 0 );
		return 1;
	}

	int oldsize = usCtgCounter;

	while ( getCntBetween ( c2, c1 ) && usCtgCounter > 1 )
	{
		usCtgCounter--;
		c2 = getTwinCtg ( upstreamCTG[usCtgCounter] );
	}

	if ( usCtgCounter != oldsize )
	{
		unsigned int prev_c2 = upstreamCTG[usCtgCounter + 1];
		unsigned int bal_prev_c2 = getTwinCtg ( prev_c2 );
		setConnectMask ( bal_prev_c2, c2, 1 );
		setConnectMask ( bal_prev_c2, c2, 0 );
		int i = usCtgCounter + 1;

		for ( ; i <= oldsize; i++ )
		{
			contig_array[upstreamCTG[i]].from_vt = prev_c2;
			contig_array[getTwinCtg ( upstreamCTG[i] )].to_vt = bal_prev_c2;
		}

		if ( ( cn_temp = getCntBetween ( c1, c2 ) ) != NULL )
		{
			setConnectMask ( c1, c2, 0 );
			setConnectDelete ( c1, c2, 0, 0 );
			return 1;
		}
	}

	len = ( min_dis + max_dis ) / 2 >= 0 ? ( min_dis + max_dis ) / 2 : 0;
	cn_temp = allocateCN ( c2, len );

	if ( cntLookupTable )
		{ putCnt2LookupTable ( c1, cn_temp ); }

	cn_temp->weight = 0;  // special connect from the original graph
	cn_temp->next = contig_array[c1].downwardConnect;
	contig_array[c1].downwardConnect = cn_temp;
	bal_c1 = getTwinCtg ( c1 );
	bal_c2 = getTwinCtg ( c2 );
	cn_temp = allocateCN ( bal_c1, len );

	if ( cntLookupTable )
		{ putCnt2LookupTable ( bal_c2, cn_temp ); }

	cn_temp->weight = 0;  // special connect from the original graph
	cn_temp->next = contig_array[bal_c2].downwardConnect;
	contig_array[bal_c2].downwardConnect = cn_temp;
	return 1;
}


/*************************************************
Function:
    catUsDsContig
Description:
    Concatenates upstream and downstream contig arrays together.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void catUsDsContig()
{
	int i;

	for ( i = 0; i < dsCtgCounter; i++ )
		{ * ( unsigned int * ) darrayPut ( solidArray, i ) = downstreamCTG[i]; }

	for ( i = usCtgCounter; i >= 0; i-- )
	{
		* ( unsigned int * ) darrayPut ( solidArray, dsCtgCounter++ ) = getTwinCtg ( upstreamCTG[i] );
	}

	solidCounter = dsCtgCounter;
}


/*************************************************
Function:
    consolidate
Description:
    Constructs scaffolds by binding connection relation among contigs
    in "solidArray".
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void consolidate()
{
	int i, j;
	CONNECT * prevCNT = NULL;
	CONNECT * cnt;
	unsigned int to_ctg;
	unsigned int from_ctg = * ( unsigned int * ) darrayGet ( solidArray, 0 );

	for ( i = 1; i < solidCounter; i++ )
	{
		to_ctg = * ( unsigned int * ) darrayGet ( solidArray, i );
		cnt = checkConnect ( from_ctg, to_ctg );

		if ( !cnt )
		{
			fprintf ( stderr, "consolidate A: no connect from %d to %d\n",
			          from_ctg, to_ctg );

			for ( j = 0; j < solidCounter; j++ )
				{ fprintf ( stderr, "%d-->", * ( unsigned int * ) darrayGet ( solidArray, j ) ); }

			fprintf ( stderr, "\n" );
			return;
		}

		cnt->singleInScaf = solidCounter == 2 ? 1 : 0;

		if ( prevCNT )
		{
			setNextInScaf ( prevCNT, cnt );
			setPrevInScaf ( cnt, 1 );
		}

		prevCNT = cnt;
		from_ctg = to_ctg;
	}

	//the reverse complementary path
	from_ctg = getTwinCtg ( * ( unsigned int * ) darrayGet ( solidArray, solidCounter - 1 ) );
	prevCNT = NULL;

	for ( i = solidCounter - 2; i >= 0; i-- )
	{
		to_ctg = getTwinCtg ( * ( unsigned int * ) darrayGet ( solidArray, i ) );
		cnt = checkConnect ( from_ctg, to_ctg );

		if ( !cnt )
		{
			fprintf ( stderr, "consolidate B: no connect from %d to %d\n", from_ctg, to_ctg );
			return;
		}

		cnt->singleInScaf = solidCounter == 2 ? 1 : 0;

		if ( prevCNT )
		{
			setNextInScaf ( prevCNT, cnt );
			setPrevInScaf ( cnt, 1 );
		}

		prevCNT = cnt;
		from_ctg = to_ctg;
	}
}

static void debugging1 ( unsigned int ctg1, unsigned int ctg2 )
{
	CONNECT * cn1;
	cn1 = getCntBetween ( ctg1, ctg2 );

	if ( cn1 )
	{
		fprintf ( stderr, "(%d,%d) mask %d deleted %d w %d,singleInScaf %d\n",
		          ctg1, ctg2, cn1->mask, cn1->deleted, cn1->weight, cn1->singleInScaf );

		if ( cn1->nextInScaf )
			{ fprintf ( stderr, "%d->%d->%d\n", ctg1, ctg2, cn1->nextInScaf->contigID ); }

		if ( cn1->prevInScaf )
			{ fprintf ( stderr, "*->%d->%d\n", ctg1, ctg2 ); }
		else if ( !cn1->nextInScaf )
			{ fprintf ( stderr, "NULL->%d->%d->NULL\n", ctg1, ctg2 ); }
	}
	else
		{ fprintf ( stderr, "%d -X- %d\n", ctg1, ctg2 ); }
}

/*************************************************
Function:
    removeTransitive
Description:
    For three contigs A, B and C, connections A->B and A->C exist.
    If there is a connection path connecting B and C or a new
    connection can be created between B and C, then deletes
    connection A->C.

       ------->
    A -> B -?-> C
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void removeTransitive()
{
	unsigned int i, bal_ctg;
	int flag = 1, out_num, in_num, count, min, max, linear;
	CONNECT * cn_temp, *cn1 = NULL, *cn2 = NULL;
	int multi_out = 0, single_out = 0, two_out = 0, may_transitive = 0, not_transitive = 0, cycle_num = 0, mask_ctg = 0, no_out = 0;
	fprintf ( stderr, "Start to remove transitive connection.\n" );

	while ( flag )
	{
		flag = 0;
		two_out = 0;
		not_transitive = 0;
		may_transitive = 0;
		cycle_num++;

		for ( i = 1; i <= num_ctg; i++ )
		{
			if ( contig_array[i].mask )
			{
				if ( cycle_num == 1 )
					{ mask_ctg++; }

				continue;
			}

			out_num = validConnect ( i, NULL );

			if ( out_num != 2 )
			{
				if ( out_num == 1 && cycle_num == 1 )
					{ single_out++; }

				if ( out_num > 2 && cycle_num == 1 )
					{ multi_out++; }

				if ( out_num == 0 && cycle_num == 1 )
					{ no_out++; }

				continue;
			}

			two_out++;
			cn_temp = contig_array[i].downwardConnect;
			count = 0;

			while ( cn_temp )
			{
				if ( cn_temp->deleted || cn_temp->mask )
				{
					cn_temp = cn_temp->next;
					continue;
				}

				count++;

				if ( count == 1 )
					{ cn1 = cn_temp; }
				else if ( count == 2 )
				{
					cn2 = cn_temp;
				}
				else
					{ break; }

				cn_temp = cn_temp->next;
			}

			if ( count > 2 )
			{
				fprintf ( stderr, "%d valid connections from ctg %d\n", count, i );
				continue;
			}

			if ( cn1->gapLen > cn2->gapLen )
			{
				cn_temp = cn1;
				cn1 = cn2;
				cn2 = cn_temp;
			}  //make sure cn1 is closer to contig i than cn2

			if ( cn1->prevInScaf && cn2->prevInScaf )
				{ continue; }

			bal_ctg = getTwinCtg ( cn2->contigID );
			in_num = validConnect ( bal_ctg, NULL );

			if ( in_num > 2 )
				{ continue; }

			int bal_c1 = getTwinCtg ( cn1->contigID );
			in_num = validConnect ( bal_c1, NULL );

			if ( in_num > 1 )
				{ continue; }

			min = cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length - ins_size_var / 2;
			max = cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length + ins_size_var / 2;

			if ( max < 0 )
				{ continue; }

			may_transitive++;
			//temprarily delete cn2
			setConnectDelete ( i, cn2->contigID, 1, 0 );
			int oldc2 = cn2->contigID;
			linear = linearC2C ( i, cn1, cn2->contigID, min, max );

			if ( linear != 1 )
			{
				not_transitive++;
				setConnectDelete ( i, cn2->contigID, 0, 0 );
				continue;
			}
			else
			{
				downstreamCTG[0] = i;
				catUsDsContig();

				if ( !checkSimple ( solidArray, solidCounter ) )
					{ continue; }

				cn1 = getCntBetween ( * ( unsigned int * ) darrayGet ( solidArray, solidCounter - 2 ),
				                      * ( unsigned int * ) darrayGet ( solidArray, solidCounter - 1 ) );

				if ( cn1 && cn1->nextInScaf && cn2->nextInScaf )
				{
					setConnectDelete ( i, cn2->contigID, 0, 0 );
					continue;
				}

				consolidate();

				if ( cn2->prevInScaf )
					substitueDSinScaf ( cn2, * ( unsigned int * ) darrayGet ( solidArray, 0 ),
					                    * ( unsigned int * ) darrayGet ( solidArray, 1 ) );

				if ( cn2->nextInScaf )
					{ substitueUSinScaf ( cn2, * ( unsigned int * ) darrayGet ( solidArray, solidCounter - 2 ) ); }

				flag++;
			}
		}

		if ( cycle_num == 1 )
		{
			fprintf ( stderr, "Total contigs                         %u\n", num_ctg );
			fprintf ( stderr, "Masked contigs                        %d\n", mask_ctg );
			fprintf ( stderr, "Remained contigs                      %u\n", num_ctg - mask_ctg );
			fprintf ( stderr, "None-outgoing-connection contigs      %d", no_out );

			if ( num_ctg - mask_ctg > 0 )
			{
				fprintf ( stderr, " (%1f%%)", ( float ) no_out / ( num_ctg - mask_ctg ) * 100 );
			}

			fprintf ( stderr, "\nSingle-outgoing-connection contigs    %d\n", single_out );
			fprintf ( stderr, "Multi-outgoing-connection contigs     %d\n", multi_out );
		}

		fprintf ( stderr, "Cycle %d\n Two-outgoing-connection contigs     %d\n Potential transitive connections    %d\n Transitive connections              %d\n", cycle_num, two_out, may_transitive, flag );

		if ( two_out > 0 )
		{
			fprintf ( stderr, " Transitive ratio                    %.1f%%\n", ( float ) flag / two_out * 100 );
		}

		if ( flag == 0 )
			{ break; }
	}
}

static void debugging2 ( unsigned int ctg )
{
	if ( ctg > num_ctg )
	{
		return;
	}

	CONNECT * cn1 = contig_array[ctg].downwardConnect;

	while ( cn1 )
	{
		if ( cn1->nextInScaf )
			{ fprintf ( stderr, "with nextInScaf %u,", cn1->nextInScaf->contigID ); }

		if ( cn1->prevInScaf )
			{ fprintf ( stderr, "with prevInScaf," ); }

		fprintf ( stderr, "%u >> %u, weight %d, gapLen %d, mask %d deleted %d, inherit %d, singleInScaf %d, bySmall %d\n",
		          ctg, cn1->contigID, cn1->weight, cn1->gapLen, cn1->mask, cn1->deleted, cn1->inherit, cn1->singleInScaf, cn1->bySmall );
		cn1 = cn1->next;
	}
}
static void debugging()
{
	//  debugging1(13298356, 13245956);
}

/*************************************************
Function:
    simplifyCnt
Description:
    Simplifys contig graph by two operations.
    1) Removes transitive connections.
    2) Picks up local contig graph and tries to line involved
        contigs.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void simplifyCnt()
{
	removeTransitive();
	debugging();
	general_linearization ( 1 );
	debugging();
}

/*************************************************
Function:
    getIndexInArray
Description:
    Gets contig's index in array.
Input:
    1. node:        contig
Output:
    None.
Return:
    Index of contig if contig existed.
    -1 otherwise.
*************************************************/
static int getIndexInArray ( unsigned int node )
{
	int index;

	for ( index = 0; index < nodeCounter; index++ )
		if ( nodesInSub[index] == node )
			{ return index; }

	return -1;
}

/*************************************************
Function:
    putNodeIntoSubgraph
Description:
    Puts contig into sub-graph.
Input:
    1. heap:            heap to store contigs
    2. distance:        current contig's distance to base contig
    3. node:            current contig
    4. index:           contig's index in array
Output:
    None.
Return:
    0 if contig already existed in sub-graph.
    1 if operation succeeded.
    -1 if index was larger than allowed maximum sub-graph size.
*************************************************/
static boolean putNodeIntoSubgraph ( FibHeap * heap, int distance, unsigned int node, int index )
{
	int pos = getIndexInArray ( node );

	if ( pos > 0 )
	{
		return 0;
	}

	if ( index >= MaxNodeInSub )
		{ return -1; }

	insertNodeIntoHeap ( heap, distance, node );
	nodesInSub[index] = node;
	nodeDistance[index] = distance;
	return 1;
}

/*************************************************
Function:
    putChainIntoSubgraph
Description:
    Puts contigs of connection chain into sub-graph.
Input:
    1. heap:            heap to store contigs
    2. distance:        current contig's distance to base contig
    3. node:            current contig
    4. index:           contig's index in array
    5. prevC:           upstream connection of current contig
Output:
    1. index:           index of next contig to add
    2. prevC:           upstream connection of last contig in connection chain
Return:
    0 if operation of putting contig into sub-graph failed.
*************************************************/
static boolean putChainIntoSubgraph ( FibHeap * heap, int distance, unsigned int node, int * index, CONNECT * prevC )
{
	unsigned int ctg = node;
	CONNECT * nextCnt;
	boolean excep, flag;
	int counter = *index;

	while ( 1 )
	{
		nextCnt = getNextContig ( ctg, prevC, &excep );

		if ( excep || !nextCnt )
		{
			*index = counter;
			return 1;
		}

		ctg = nextCnt->contigID;
		distance += nextCnt->gapLen + contig_array[ctg].length;
		flag = putNodeIntoSubgraph ( heap, distance, ctg, counter );

		if ( flag < 0 )
			{ return 0; }

		if ( flag > 0 )
			{ counter++; }

		prevC = nextCnt;
	}
}

/*************************************************
Function:
    checkUnique
Description:
    Checks if a contig is unique by trying to line its upstream and
    downstream contigs.
Input:
    1. node:            contig to check
    2. tolerance:       allowed percentage that overlap length among
                    contigs in sub-graph accounts for in all involved contigs' length
Output:
    None.
Return:
    1 if contig was unique.
*************************************************/
static boolean checkUnique ( unsigned int node, double tolerance )
{
	CONNECT * ite_cnt;
	unsigned int currNode;
	int distance;
	int popCounter = 0;
	boolean flag;
	currNode = node;
	FibHeap * heap = newFibHeap();
	putNodeIntoSubgraph ( heap, 0, currNode, 0 );
	nodeCounter = 1;
	ite_cnt = contig_array[currNode].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->deleted || ite_cnt->mask )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		currNode = ite_cnt->contigID;
		distance = ite_cnt->gapLen + contig_array[currNode].length;
		flag = putNodeIntoSubgraph ( heap, distance, currNode, nodeCounter );

		if ( flag < 0 )
		{
			destroyHeap ( heap );
			return 0;
		}

		if ( flag > 0 )
			{ nodeCounter++; }

		flag = putChainIntoSubgraph ( heap, distance, currNode, &nodeCounter, ite_cnt );

		if ( !flag )
		{
			destroyHeap ( heap );
			return 0;
		}

		ite_cnt = ite_cnt->next;
	}

	if ( nodeCounter <= 2 ) // no more than 2 valid connections
	{
		destroyHeap ( heap );
		return 1;
	}

	while ( ( currNode = removeNextNodeFromHeap ( heap ) ) != 0 )
		{ nodesInSubInOrder[popCounter++] = currNode; }

	destroyHeap ( heap );
	flag = checkOverlapInBetween ( tolerance );
	return flag;
}

/*************************************************
Function:
    maskRepeat
Description:
    Masks repeat contigs.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void maskRepeat()
{
	int in_num, out_num, flagA, flagB;
	int counter = 0;
	int puzzleCounter = 0;
	unsigned int i, bal_i;

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].mask )
			{ continue; }

		bal_i = getTwinCtg ( i );
		in_num = validConnect ( bal_i, NULL );
		out_num = validConnect ( i, NULL );

		if ( in_num > 1 || out_num > 1 )
			{ puzzleCounter++; }
		else
		{
			if ( isSmallerThanTwin ( i ) )
				{ i++; }

			continue;
		}

		if ( contig_array[i].cvg > 1.4 * cvgAvg )
		{
			counter++;
			maskContig ( i, 1 );

			if ( isSmallerThanTwin ( i ) )
				{ i++; }

			continue;
		}

		if ( in_num > 1 )
			{ flagA = checkUnique ( bal_i, OverlapPercent ); }
		else
			{ flagA = 1; }

		if ( out_num > 1 )
			{ flagB = checkUnique ( i, OverlapPercent ); }
		else
			{ flagB = 1; }

		if ( !flagA || !flagB )
		{
			counter++;
			maskContig ( i, 1 );
		}

		if ( isSmallerThanTwin ( i ) )
			{ i++; }
	}

	fprintf ( stderr, "Mask repeats:\n Puzzles           %d\n Masked contigs    %d\n", puzzleCounter, counter );
}

/*************************************************
Function:
    Countlink
Description:
    Counts valid connections of all contigs.
Input:
    None.
Output:
    None.
Return:
    Valid connections of all contigs.
*************************************************/
static int Countlink()
{
	unsigned int i, bal_i;
	int conflict_count = 0;

	for ( i = 1; i < num_ctg; i++ )
	{
		if ( contig_array[i].mask )
			{ continue; }

		int out_num = validConnect ( i, NULL );

		if ( out_num > 1 )
			{ conflict_count++; }
	}

	return conflict_count;
}

/*************************************************
Function:
    ordering
Description:
    Constructs scaffold using alignment information of paired-end
    reads rank by rank.
Input:
    1. deWeak:      indicator of deleting weak connection.
    2. downS:       indicator of large insert size paired-end reads
    3. nonlinear:       indicator of simplifying graph using loose restraint
    4. infile:          prefix of graph
Output:
    None.
Return:
    None.
*************************************************/
static void ordering ( boolean deWeak, boolean downS, boolean nonlinear, char * infile )
{
	int conf0, conf1, conf2, conf3, conf4, conf5;
	debugging();

	if ( downS )
	{
		downSlide();
		debugging();

		if ( deWeak )
			{ deleteWeakCnt ( weakPE ); }
	}
	else
	{
		if ( deWeak )
			{ deleteWeakCnt ( weakPE ); }
	}

	debugging();
	simplifyCnt();
	debugging();
	maskRepeat();
	debugging();
	simplifyCnt();
	debugging();

	if ( nonlinear )
	{
		fprintf ( stderr, "Non-strict linearization.\n" );
		general_linearization ( 0 );
	}

	maskPuzzle ( 2, 0 );
	debugging();
	freezing();
	debugging();
}


/*************************************************
Function:
    checkOverlapInBetween
Description:
    Checks if adjacent contigs in the array have reasonable overlap.
Input:
    1. tolerance:       max percentage that overlap length accounts for
Output:
    None.
Return:
    1 if the overlap situation was resonable.
*************************************************/
boolean checkOverlapInBetween ( double tolerance )
{
	int i, gap;
	int index;
	unsigned int node;
	int lenSum, lenOlp;
	lenSum = lenOlp = 0;

	for ( i = 0; i < nodeCounter; i++ )
	{
		node = nodesInSubInOrder[i];
		lenSum += contig_array[node].length;
		index = getIndexInArray ( node );
		nodeDistanceInOrder[i] = nodeDistance[index];
	}

	if ( lenSum < 1 )
		{ return 1; }

	for ( i = 0; i < nodeCounter - 1; i++ )
	{
		gap = nodeDistanceInOrder[i + 1] - nodeDistanceInOrder[i]
		      - contig_array[nodesInSubInOrder[i + 1]].length;

		if ( -gap > 0 )
			{ lenOlp += -gap; }

		if ( ( double ) lenOlp / lenSum > tolerance )
		{
			return 0;
		}
	}

	return 1;
}


/*********   the following codes are for freezing current scaffolds   ****************/
/*************************************************
Function:
    setUsed
Description:
    Checks status of connection between contigs. If none of them
    equals to "flag", sets all status to "flag" and adjusts the related
    indicators accordingly. Otherwise does NOT change the status
    at all.
Input:
    1. start:           the first contig
    2. array:           contig array
    3. max_steps:       max number of connection to check
    4. flag:            new status
Output:
    None.
Return:
    o if setting successed.
*************************************************/
static boolean setUsed ( unsigned int start, unsigned int * array, int max_steps, boolean flag )
{
	unsigned int prevCtg = start;
	unsigned int twinA, twinB;
	int j;
	CONNECT * cnt;
	boolean usedFlag = 0;
	// save 'used' to 'checking'
	prevCtg = start;

	for ( j = 0; j < max_steps; j++ )
	{
		if ( array[j] == 0 )
			{ break; }

		cnt = getCntBetween ( prevCtg, array[j] );

		if ( !cnt )
		{
			fprintf ( stderr, "setUsed: no connect between %d and %d\n", prevCtg, array[j] );
			prevCtg = array[j];
			continue;
		}

		if ( cnt->used == flag || cnt->nextInScaf || cnt->prevInScaf || cnt->singleInScaf )
		{
			return 1;
		}

		cnt->checking = cnt->used;
		twinA = getTwinCtg ( prevCtg );
		twinB = getTwinCtg ( array[j] );
		cnt = getCntBetween ( twinB, twinA );

		if ( cnt )
			{ cnt->checking = cnt->used; }

		prevCtg = array[j];
	}

	// set used to flag
	prevCtg = start;

	for ( j = 0; j < max_steps; j++ )
	{
		if ( array[j] == 0 )
			{ break; }

		cnt = getCntBetween ( prevCtg, array[j] );

		if ( !cnt )
		{
			prevCtg = array[j];
			continue;
		}

		if ( cnt->used == flag )
		{
			usedFlag = 1;
			break;
		}

		cnt->used = flag;
		twinA = getTwinCtg ( prevCtg );
		twinB = getTwinCtg ( array[j] );
		cnt = getCntBetween ( twinB, twinA );

		if ( cnt )
			{ cnt->used = flag; }

		prevCtg = array[j];
	}

	// set mask to 'NOT flag' or set used to original value
	prevCtg = start;

	for ( j = 0; j < max_steps; j++ )
	{
		if ( array[j] == 0 )
			{ break; }

		cnt = getCntBetween ( prevCtg, array[j] );

		if ( !cnt )
		{
			prevCtg = array[j];
			continue;
		}

		if ( !usedFlag )
			{ cnt->mask = 1 - flag; }
		else
			{ cnt->used = cnt->checking; }

		twinA = getTwinCtg ( prevCtg );
		twinB = getTwinCtg ( array[j] );
		cnt = getCntBetween ( twinB, twinA );
		cnt->used = 1 - flag;

		if ( !usedFlag )
			{ cnt->mask = 1 - flag; }
		else
			{ cnt->used = cnt->checking; }

		prevCtg = array[j];
	}

	return usedFlag;
}


/*************************************************
Function:
    score_pass
Description:
    Checks whether a contig is supposed to locate in a scaffold
    according to its connections to other contigs in this scaffold.
Input:
    1. array:           array of contigs of part of or whole scaffold
    2. Counter:     contig number in array
    3. beforep:     index of previous contig
    4. afterp:          index of next contig
    5. id:          current contig
Output:
    None.
Return:
    1 if this contig had enough support from other contigs.
*************************************************/
int score_pass ( DARRAY * array, int Counter, int beforep, int afterp, int id )
{
	int outnum = allConnect ( id, NULL );
	int innum = allConnect ( getTwinCtg ( id ), NULL );
	int start = beforep - 2 * innum > 0 ? beforep - 2 * innum : 0;
	int end = afterp + 2 * outnum < Counter ? afterp + 2 * outnum : Counter;
	int i, inc = 1, outc = 1;
	CONNECT * dh_cnt;

	for ( i = start; i < end; i++ )
	{
		if ( i < beforep )
		{
			dh_cnt = getCntBetween ( * ( unsigned int * ) darrayGet ( array, i ), id );

			if ( dh_cnt )
				{ inc++; }
		}

		if ( i > afterp )
		{
			dh_cnt = getCntBetween ( id, * ( unsigned int * ) darrayGet ( array, i ) );

			if ( dh_cnt )
				{ outc++; }
		}
	}

	if ( inc == innum || outc == outnum )
		{ return 1; }

	if ( ( inc == 1 && innum > 2 ) || ( outc == 1 && outnum > 2 ) )
		{ return 0; }

	int score = ( int ) ( ( ( double ) ( inc * outc ) / ( double ) ( innum * outnum ) ) * 100 );

	if ( score > 30 )
		{ return 1; }

	return 0;
}

/*************************************************
Function:
    recoverMask
Description:
    Searchs a contig path between two contigs accroding to
    connection relation. If a path is found, recovers those masked
    connection between contigs in this path.
    A  ------  D
       > [i] <
    B          E
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void recoverMask()
{
	unsigned int i, ctg, bal_ctg, start, finish;
	int num3, num5, j, t;
	CONNECT * bindCnt, *cnt;
	int min, max, max_steps = 5, num_route, length;
	int tempCounter, recoverCounter = 0;
	boolean multiUSE, change;
	int stat[] = {0, 0, 0, 0, 0, 0, 0};

	for ( i = 1; i <= num_ctg; i++ )
		{ contig_array[i].flag = 0; }

	so_far = ( unsigned int * ) ckalloc ( max_n_routes * sizeof ( unsigned int ) );
	found_routes = ( unsigned int ** ) ckalloc ( max_n_routes * sizeof ( unsigned int * ) );

	for ( j = 0; j < max_n_routes; j++ )
		{ found_routes[j] = ( unsigned int * ) ckalloc ( max_steps * sizeof ( unsigned int ) ); }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].flag || contig_array[i].mask || !contig_array[i].downwardConnect )
			{ continue; }

		bindCnt = getBindCnt ( i );

		if ( !bindCnt )
			{ continue; }

		//first scan get the average coverage by longer pe
		num5 = num3 = 0;
		ctg = i;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg ( i )].flag = 1;

		while ( bindCnt )
		{
			if ( bindCnt->used )
				{ break; }

			setConnectUsed ( ctg, bindCnt->contigID, 1 );
			ctg = bindCnt->contigID;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg ( i );
		bindCnt = getBindCnt ( ctg );

		while ( bindCnt )
		{
			if ( bindCnt->used )
				{ break; }

			setConnectUsed ( ctg, bindCnt->contigID, 1 );
			ctg = bindCnt->contigID;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}

		if ( num5 + num3 < 2 )
			{ continue; }

		tempCounter = solidCounter = 0;

		for ( j = num3 - 1; j >= 0; j-- )
			* ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf3, j );

		for ( j = 0; j < num5; j++ )
			* ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf5, j );

		change = 0;

		for ( t = 0; t < tempCounter - 1; t++ )
		{
			* ( unsigned int * ) darrayPut ( solidArray, solidCounter++ ) =
			    * ( unsigned int * ) darrayGet ( tempArray, t );
			start = * ( unsigned int * ) darrayGet ( tempArray, t );
			finish = * ( unsigned int * ) darrayGet ( tempArray, t + 1 );
			num_route = num_trace = 0;
			cnt = checkConnect ( start, finish );

			if ( !cnt )
			{
				fprintf ( stderr, "Warning from recoverMask: no connection (%d %d), start at %d\n",
				          start, finish, i );
				cnt = getCntBetween ( start, finish );

				if ( cnt )
					{ debugging1 ( start, finish ); }

				continue;
			}

			length = cnt->gapLen + contig_array[finish].length;
			min = length - 1.5 * ins_size_var;
			max = length + 1.5 * ins_size_var;
			traceAlongMaskedCnt ( finish, start, max_steps, min, max, 0, 0, &num_route );

			if ( finish == start )
			{
				for ( j = 0; j < tempCounter; j++ )
					{ fprintf ( stderr, "->%d", * ( unsigned int * ) darrayGet ( tempArray, j ) ); }

				fprintf ( stderr, ": start at %d\n", i );
			}

			if ( num_route == 1 )
			{
				for ( j = 0; j < max_steps; j++ )
					if ( found_routes[0][j] == 0 )
						{ break; }

				if ( j < 1 )
					{ continue; }

				//check if connects have been used more than once
				multiUSE = setUsed ( start, found_routes[0], max_steps, 1 );

				if ( multiUSE )
					{ continue; }

				for ( j = 0; j < max_steps; j++ )
				{
					if ( j + 1 == max_steps || found_routes[0][j + 1] == 0 )
						{ break; }

					* ( unsigned int * ) darrayPut ( solidArray, solidCounter++ ) = found_routes[0][j];
					contig_array[found_routes[0][j]].flag = 1;
					contig_array[getTwinCtg ( found_routes[0][j] )].flag = 1;
				}

				recoverCounter += j;
				setConnectDelete ( start, finish, 1, 1 );
				change = 1;
				stat[0]++;
			}  //end if num_route=1
			else if ( num_route > 1 )
			{
				//multi-route.
				int k, l, num, longest = 0, longestid = 0;
				int merg = 0, quality = 0;

				// get the longest route.
				for ( k = 0; k < num_route; k++ )
				{
					for ( j = 0; j < max_steps; j++ )
					{
						if ( j + 1 == max_steps || found_routes[k][j + 1] == 0 )
						{
							if ( j > longest )
							{
								longest = j;
								longestid = k;
							}

							break;
						}
					}
				}

				stat[1]++;

				if ( longest == 1 ) //multi one.
				{
					stat[2]++;

					for ( k = 0; k < num_route; k++ )
					{
						if ( score_pass ( tempArray, tempCounter, t, t + 1, found_routes[k][0] ) )
						{
							longestid = k;
							quality = 1;
							stat[3]++;
							break;
						}
					}

					if ( quality == 0 )
					{
						continue;
					}
				}
				else
				{
					stat[4]++;

					for ( k = 0; k < num_route; k++ )
					{
						if ( k == longestid )
							{ continue; }

						int merg_num = 0, total = 0;

						for ( j = 0; j < max_steps; j++ )
						{
							if ( j + 1 == max_steps || found_routes[k][j + 1] == 0 )
							{
								total = j;
								break;
							}

							for ( l = 0; l < longest; l++ )
							{
								if ( found_routes[k][j] == found_routes[longestid][l] )
								{
									merg_num++;
									break;
								}
							}
						}

						if ( merg_num == total )
							{ merg++; }
					}
				}

				if ( merg == num_route - 1 || quality == 1 || merg >= longest )
				{
					multiUSE = setUsed ( start, found_routes[longestid], max_steps, 1 );

					if ( multiUSE )
						{ continue; }

					stat[5]++;

					for ( j = 0; j < longest; j++ )
					{
						* ( unsigned int * ) darrayPut ( solidArray, solidCounter++ ) = found_routes[longestid][j];
						contig_array[found_routes[longestid][j]].flag = 1;
						contig_array[getTwinCtg ( found_routes[longestid][j] )].flag = 1;
					}

					stat[6] += j;
					recoverCounter += j;
					setConnectDelete ( start, finish, 1, 1 );
					change = 1;
				}
			}
		}

		* ( unsigned int * ) darrayPut ( solidArray, solidCounter++ ) =
		    * ( unsigned int * ) darrayGet ( tempArray, tempCounter - 1 );

		if ( change )
			{ consolidate(); }
	}

	fprintf ( stderr, "\nRecover contigs.\n" );
	fprintf ( stderr, " Total recovered contigs    %d\n", recoverCounter );
	fprintf ( stderr, " Single-route cases         %d\n", stat[0] );
	fprintf ( stderr, " Multi-route cases          %d\n", stat[1] );

	for ( i = 1; i <= num_ctg; i++ )
	{
		cnt = contig_array[i].downwardConnect;

		while ( cnt )
		{
			cnt->used = 0;
			cnt->checking = 0;
			cnt = cnt->next;
		}
	}

	for ( j = 0; j < max_n_routes; j++ )
		{ free ( ( void * ) found_routes[j] ); }

	free ( ( void * ) found_routes );
	free ( ( void * ) so_far );
}


/*************************************************
Function:
    unBindLink
Description:
    Destroys upstream and downstram connections of connection
    between two specified contigs.
Input:
    1. CB:      the first specified contig
    2. CC:      the second specified contig
Output:
    None.
Return:
    None.
*************************************************/
static void unBindLink ( unsigned int CB, unsigned int CC )
{
	CONNECT * cnt1 = getCntBetween ( CB, CC );

	if ( !cnt1 )
		{ return; }

	if ( cnt1->singleInScaf )
		{ cnt1->singleInScaf = 0; }

	CONNECT * cnt2 = getCntBetween ( getTwinCtg ( CC ), getTwinCtg ( CB ) );

	if ( !cnt2 )
		{ return; }

	if ( cnt2->singleInScaf )
		{ cnt2->singleInScaf = 0; }

	if ( cnt1->nextInScaf )
	{
		unsigned int CD = cnt1->nextInScaf->contigID;
		cnt1->nextInScaf->prevInScaf = 0;
		cnt1->nextInScaf = NULL;
		CONNECT * cnt3 = getCntBetween ( getTwinCtg ( CD ), getTwinCtg ( CC ) );

		if ( cnt3 )
			{ cnt3->nextInScaf = NULL; }

		cnt2->prevInScaf = 0;
	}

	if ( cnt2->nextInScaf )
	{
		unsigned int bal_CA = cnt2->nextInScaf->contigID;
		cnt2->nextInScaf->prevInScaf = 0;
		cnt2->nextInScaf = NULL;
		CONNECT * cnt4 = getCntBetween ( getTwinCtg ( bal_CA ), CB );

		if ( cnt4 )
			{ cnt4->nextInScaf = NULL; }

		cnt1->prevInScaf = 0;
	}
}

/*************************************************
Function:
    freezing
Description:
    Updates connection relation and scaffold ID of contigs in a scaffold.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void freezing()
{
	int num5, num3;
	unsigned int ctg, bal_ctg;
	unsigned int i;
	int j, t;
	CONNECT * cnt, *prevCNT, *nextCnt;
	boolean excep;

	for ( i = 1; i <= num_ctg; i++ )
	{
		contig_array[i].flag = 0;
		contig_array[i].from_vt = 0;
		contig_array[i].to_vt = 0;
		cnt = contig_array[i].downwardConnect;

		while ( cnt )
		{
			cnt->used = 0;
			cnt->checking = 0;
			cnt->singleInScaf = 0;
			cnt = cnt->next;
		}
	}

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].flag || contig_array[i].mask )
			{ continue; }

		if ( !contig_array[i].downwardConnect || !validConnect ( i, NULL ) )
		{
			continue;
		}

		num5 = num3 = 0;
		ctg = i;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg ( i )].flag = 1;
		prevCNT = NULL;
		cnt = getNextContig ( ctg, prevCNT, &excep );

		while ( cnt )
		{
			if ( contig_array[cnt->contigID].flag )
			{
				unBindLink ( ctg, cnt->contigID );
				break;
			}

			nextCnt = getNextContig ( cnt->contigID, cnt, &excep );
			setConnectUsed ( ctg, cnt->contigID, 1 );
			ctg = cnt->contigID;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		ctg = getTwinCtg ( i );

		if ( num5 >= 2 )
			{ prevCNT = checkConnect ( getTwinCtg ( * ( unsigned int * ) darrayGet ( scaf5, 1 ) ), ctg ); }
		else
			{ prevCNT = NULL; }

		cnt = getNextContig ( ctg, prevCNT, &excep );

		while ( cnt )
		{
			if ( contig_array[cnt->contigID].flag )
			{
				unBindLink ( ctg, cnt->contigID );
				break;
			}

			nextCnt = getNextContig ( cnt->contigID, cnt, &excep );
			setConnectUsed ( ctg, cnt->contigID, 1 );
			ctg = cnt->contigID;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = bal_ctg;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		if ( num5 + num3 < 2 )
			{ continue; }

		solidCounter = 0;

		for ( j = num3 - 1; j >= 0; j-- )
			* ( unsigned int * ) darrayPut ( solidArray, solidCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf3, j );

		for ( j = 0; j < num5; j++ )
			* ( unsigned int * ) darrayPut ( solidArray, solidCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf5, j );

		unsigned int firstCtg = 0;
		unsigned int lastCtg = 0;
		unsigned int firstTwin = 0;
		unsigned int lastTwin = 0;

		for ( t = 0; t < solidCounter; t++ )
			if ( !contig_array[* ( unsigned int * ) darrayGet ( solidArray, t )].mask )
			{
				firstCtg = * ( unsigned int * ) darrayGet ( solidArray, t );
				break;
			}

		for ( t = solidCounter - 1; t >= 0; t-- )
			if ( !contig_array[* ( unsigned int * ) darrayGet ( solidArray, t )].mask )
			{
				lastCtg = * ( unsigned int * ) darrayGet ( solidArray, t );
				break;
			}

		if ( firstCtg == 0 || lastCtg == 0 )
		{
			fprintf ( stderr, "scaffold start at %d, stop at %d, freezing began with %d\n", firstCtg, lastCtg, i );

			for ( j = 0; j < solidCounter; j++ )
				fprintf ( stderr, "->%d(%d %d)", * ( unsigned int * ) darrayGet ( solidArray, j )
				          , contig_array[* ( unsigned int * ) darrayGet ( solidArray, j )].mask
				          , contig_array[* ( unsigned int * ) darrayGet ( solidArray, j )].flag );

			fprintf ( stderr, "\n" );
		}
		else
		{
			firstTwin = getTwinCtg ( firstCtg );
			lastTwin = getTwinCtg ( lastCtg );
		}

		for ( t = 0; t < solidCounter; t++ )
		{
			unsigned int ctg = * ( unsigned int * ) darrayGet ( solidArray, t );

			if ( contig_array[ctg].from_vt > 0 )
			{
				contig_array[ctg].mask = 1;
				contig_array[getTwinCtg ( ctg )].mask = 1;
				fprintf ( stderr, "Repeat: contig %d (%d) appears more than once\n", ctg, getTwinCtg ( ctg ) );
			}
			else
			{
				contig_array[ctg].from_vt = firstCtg;
				contig_array[ctg].to_vt = lastCtg;
				contig_array[ctg].indexInScaf = t + 1;
				contig_array[getTwinCtg ( ctg )].from_vt = lastTwin;
				contig_array[getTwinCtg ( ctg )].to_vt = firstTwin;
				contig_array[getTwinCtg ( ctg )].indexInScaf = solidCounter - t;
			}
		}

		consolidate();
	}

	fprintf ( stderr, "Freezing done.\n" );
	fflush ( stdout );

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].flag )
			{ contig_array[i].flag = 0; }

		if ( contig_array[i].from_vt == 0 )
		{
			contig_array[i].from_vt = i;
			contig_array[i].to_vt = i;
		}

		cnt = contig_array[i].downwardConnect;

		while ( cnt )
		{
			cnt->used = 0;
			cnt->checking = 0;
			cnt = cnt->next;
		}
	}
}

/************** codes below this line are for pulling the scaffolds out ************/
/*************************************************
Function:
    output1gap
Description:
    Outputs information of a filled gap.
Input:
    1. fo:          output file
    2. max_steps:       maximum allowed found contigs for filling gap
Output:
    None.
Return:
    None.
*************************************************/
void output1gap ( FILE * fo, int max_steps )
{
	int i, len, seg;
	len = seg = 0;

	for ( i = 0; i < max_steps - 1; i++ )
	{
		if ( found_routes[0][i + 1] == 0 )
			{ break; }

		len += contig_array[found_routes[0][i]].length;
		seg++;
	}

	fprintf ( fo, "GAP %d %d", len, seg );

	for ( i = 0; i < max_steps - 1; i++ )
	{
		if ( found_routes[0][i + 1] == 0 )
			{ break; }

		fprintf ( fo, " %d", found_routes[0][i] );
	}

	fprintf ( fo, "\n" );
}

static int weakCounter;

/*************************************************
Function:
    printCnts
Description:
    Outputs connection information of specified contig.
Input:
    1. fp:      output file
    2. ctg:     specified contig
Output:
    None.
Return:
    0 if contig's downstream connection was not weak connection.
*************************************************/
static boolean printCnts ( FILE * fp, unsigned int ctg )
{
	CONNECT * cnt = contig_array[ctg].downwardConnect;
	boolean flag = 0, ret = 0;
	unsigned int bal_ctg = getTwinCtg ( ctg );
	unsigned int linkCtg;

	if ( isSameAsTwin ( ctg ) )
		{ return ret; }

	CONNECT * bindCnt = getBindCnt ( ctg );

	if ( bindCnt && bindCnt->bySmall && bindCnt->weakPoint )
	{
		weakCounter++;
		fprintf ( fp, "\tWP" );
		ret = 1;
	}

	while ( cnt )
	{
		if ( cnt->weight && !cnt->inherit )
		{
			if ( !flag )
			{
				flag = 1;
				fprintf ( fp, "\t#DOWN " );
			}

			linkCtg = cnt->contigID;

			if ( isLargerThanTwin ( linkCtg ) )
				{ linkCtg = getTwinCtg ( linkCtg ); }

			fprintf ( fp, "%d:%d:%d ", index_array[linkCtg], cnt->weight, cnt->gapLen );
		}

		cnt = cnt->next;
	}

	flag = 0;
	cnt = contig_array[bal_ctg].downwardConnect;

	while ( cnt )
	{
		if ( cnt->weight && !cnt->inherit )
		{
			if ( !flag )
			{
				flag = 1;
				fprintf ( fp, "\t#UP " );
			}

			linkCtg = cnt->contigID;

			if ( isLargerThanTwin ( linkCtg ) )
				{ linkCtg = getTwinCtg ( linkCtg ); }

			fprintf ( fp, "%d:%d:%d ", index_array[linkCtg], cnt->weight, cnt->gapLen );
		}

		cnt = cnt->next;
	}

	fprintf ( fp, "\n" );
	return ret;
}

/*************************************************
Function:
    ScafStat
Description:
    Makes statistic of scaffolds and contigs, including total length,
    non-N length, number, N50, etc.
Input:
    1. len_cut:     length cutoff
    2. graphfile:       prefix of graph
Output:
    None.
Return:
    None.
*************************************************/
void ScafStat ( int len_cut, char * graphfile )
{
	FILE * fp, *fp2, *fo;
	char line[1024];
	sprintf ( line, "%s.scafSeq", graphfile );
	fp = ckopen ( line, "r" );
	sprintf ( line, "%s.contig", graphfile );
	fp2 = ckopen ( line, "r" );
	sprintf ( line, "%s.scafStatistics", graphfile );
	fo = ckopen ( line, "w" );
	fprintf ( fo, "<-- Information for assembly Scaffold '%s.scafSeq'.(cut_off_length < 100bp) -->\n\n", graphfile );
	int  cut_off_len = 0;
	char Nucleotide;
	char buf[4000];
	long Scaffold_Number = 0;
	long Scaffold_Number_Scaf = 0;
	long Scaffold_Number_Contig = 0;
	long Singleton_Number_Scaf = 0;
	long Singleton_Number = 0;
	long * Singleton_Seq  = ( long * ) malloc ( sizeof ( long ) );
	long long A_num_all = 0;
	long * A_num  = ( long * ) malloc ( sizeof ( long ) );
	long long C_num_all = 0;
	long * C_num  = ( long * ) malloc ( sizeof ( long ) );
	long long G_num_all = 0;
	long * G_num  = ( long * ) malloc ( sizeof ( long ) );
	long long T_num_all = 0;
	long * T_num  = ( long * ) malloc ( sizeof ( long ) );
	long long N_num_all = 0;
	long * N_num  = ( long * ) malloc ( sizeof ( long ) );
	long long Non_ACGTN_all = 0;
	long * Non_ACGTN  = ( long * ) malloc ( sizeof ( long ) );
	long long Size_includeN = 0;
	long * Size_Seq  = ( long * ) malloc ( sizeof ( long ) );
	int  k;
	long long Sum = 0;
	int  flag[10];
	int  flag_known = 0;
	long n100 = 0;
	long n500 = 0;
	long n1k  = 0;
	long n10k = 0;
	long n100k = 0;
	long n1m  = 0;
	long N50 = 0;
	long N50_known = 0;
	long Num_N50_known = 0;
	cut_off_len = len_cut;
	A_num[Scaffold_Number] = 0;
	C_num[Scaffold_Number] = 0;
	G_num[Scaffold_Number] = 0;
	T_num[Scaffold_Number] = 0;
	N_num[Scaffold_Number] = 0;
	Non_ACGTN[Scaffold_Number] = 0;
	Size_Seq[Scaffold_Number] = 0;
	Singleton_Seq[Scaffold_Number] = 0;
	Nucleotide = fgetc ( fp );

	while ( Nucleotide != EOF )
	{
		if ( Nucleotide == '>' )
		{
			if ( ( Scaffold_Number > 0 ) && ( Size_Seq[Scaffold_Number - 1] < cut_off_len ) )
			{
				A_num_all = A_num_all - A_num[Scaffold_Number - 1];
				C_num_all = C_num_all - C_num[Scaffold_Number - 1];
				G_num_all = G_num_all - G_num[Scaffold_Number - 1];
				T_num_all = T_num_all - T_num[Scaffold_Number - 1];
				N_num_all = N_num_all - N_num[Scaffold_Number - 1];
				Non_ACGTN_all = Non_ACGTN_all - Non_ACGTN[Scaffold_Number - 1];
				Size_includeN = Size_includeN - Size_Seq[Scaffold_Number - 1];
				Singleton_Number = Singleton_Number - Singleton_Seq[Scaffold_Number - 1];
				Scaffold_Number = Scaffold_Number - 1;
			}
			else
			{
				Size_Seq = ( long * ) realloc ( Size_Seq, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				A_num    = ( long * ) realloc ( A_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				C_num    = ( long * ) realloc ( C_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				G_num    = ( long * ) realloc ( G_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				T_num    = ( long * ) realloc ( T_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				N_num    = ( long * ) realloc ( N_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				Non_ACGTN = ( long * ) realloc ( Non_ACGTN, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				Singleton_Seq = ( long * ) realloc ( Singleton_Seq, ( Scaffold_Number + 2 ) * sizeof ( long ) );
			}

			Scaffold_Number++;
			A_num[Scaffold_Number - 1] = 0;
			C_num[Scaffold_Number - 1] = 0;
			G_num[Scaffold_Number - 1] = 0;
			T_num[Scaffold_Number - 1] = 0;
			N_num[Scaffold_Number - 1] = 0;
			Non_ACGTN[Scaffold_Number - 1] = 0;
			Size_Seq[Scaffold_Number - 1] = 0;
			Singleton_Seq[Scaffold_Number - 1] = 0;
			Nucleotide = fgetc ( fp );

			if ( Nucleotide == 'C' )
			{
				Singleton_Number++;
				Singleton_Seq[Scaffold_Number - 1] ++;
			}

			fgets ( buf, 4000, fp );
		}
		else if ( ( Nucleotide == 'N' ) || ( Nucleotide == 'n' ) )
		{
			N_num[Scaffold_Number - 1] ++;
			N_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'A' ) || ( Nucleotide == 'a' ) )
		{
			A_num[Scaffold_Number - 1] ++;
			A_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'C' ) || ( Nucleotide == 'c' ) )
		{
			C_num[Scaffold_Number - 1] ++;
			C_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'G' ) || ( Nucleotide == 'g' ) )
		{
			G_num[Scaffold_Number - 1] ++;
			G_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'T' ) || ( Nucleotide == 't' ) )
		{
			T_num[Scaffold_Number - 1] ++;
			T_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else
		{
			if ( ( Nucleotide != '\n' ) && ( Nucleotide != '\r' ) )
			{
				Non_ACGTN[Scaffold_Number - 1] ++;
				Non_ACGTN_all++;
				Size_Seq[Scaffold_Number - 1] ++;
				Size_includeN++;
			}
		}

		Nucleotide = fgetc ( fp );
	}

	if ( Size_Seq[Scaffold_Number - 1] < cut_off_len )
	{
		A_num_all = A_num_all - A_num[Scaffold_Number - 1];
		C_num_all = C_num_all - C_num[Scaffold_Number - 1];
		G_num_all = G_num_all - G_num[Scaffold_Number - 1];
		T_num_all = T_num_all - T_num[Scaffold_Number - 1];
		N_num_all = N_num_all - N_num[Scaffold_Number - 1];
		Non_ACGTN_all = Non_ACGTN_all - Non_ACGTN[Scaffold_Number - 1];
		Size_includeN = Size_includeN - Size_Seq[Scaffold_Number - 1];
		Singleton_Number = Singleton_Number - Singleton_Seq[Scaffold_Number - 1];
		Scaffold_Number = Scaffold_Number - 1;
	}

	qsort ( Size_Seq, Scaffold_Number, sizeof ( Size_Seq[0] ), cmp_int );
	fprintf ( fo, "Size_includeN\t%lld\n", Size_includeN );
	fprintf ( fo, "Size_withoutN\t%lld\n", Size_includeN - N_num_all );
	fprintf ( fo, "Scaffold_Num\t%ld\n", Scaffold_Number );
	fprintf ( fo, "Mean_Size\t%lld\n", Size_includeN / Scaffold_Number );
	fprintf ( fo, "Median_Size\t%ld\n", Size_Seq[ ( Scaffold_Number + 1 ) / 2 - 1] );
	fprintf ( fo, "Longest_Seq\t%ld\n", Size_Seq[Scaffold_Number - 1] );
	fprintf ( fo, "Shortest_Seq\t%ld\n", Size_Seq[0] );
	fprintf ( fo, "Singleton_Num\t%ld\n", Singleton_Number );
	fprintf ( fo, "Average_length_of_break(N)_in_scaffold\t%lld\n", N_num_all / Scaffold_Number );
	fprintf ( fo, "\n" );

	if ( known_genome_size )
	{
		fprintf ( fo, "Known_genome_size\t%ld\n", known_genome_size );
		fprintf ( fo, "Total_scaffold_length_as_percentage_of_known_genome_size\t%.2f%\n", 100 * ( 1.0 * Size_includeN / known_genome_size ) );
	}
	else
	{
		fprintf ( fo, "Known_genome_size\tNaN\n" );
		fprintf ( fo, "Total_scaffold_length_as_percentage_of_known_genome_size\tNaN\n" );
	}

	fprintf ( fo, "\n" );

	for ( k = 0; k < Scaffold_Number; k++ )
	{
		if ( Size_Seq[k] > 100 )
		{
			n100++;
		}

		if ( Size_Seq[k] > 500 )
		{
			n500++;
		}

		if ( Size_Seq[k] > 1000 )
		{
			n1k++;
		}

		if ( Size_Seq[k] > 10000 )
		{
			n10k++;
		}

		if ( Size_Seq[k] > 100000 )
		{
			n100k++;
		}

		if ( Size_Seq[k] > 1000000 )
		{
			n1m++;
		}
	}

	fprintf ( fo, "scaffolds>100 \t%ld\t%.2f%\n", n100 , 100 * ( 1.0 * n100 / Scaffold_Number ) );
	fprintf ( fo, "scaffolds>500 \t%ld\t%.2f%\n", n500 , 100 * ( 1.0 * n500 / Scaffold_Number ) );
	fprintf ( fo, "scaffolds>1K  \t%ld\t%.2f%\n", n1k   , 100 * ( 1.0 * n1k / Scaffold_Number ) );
	fprintf ( fo, "scaffolds>10K \t%ld\t%.2f%\n", n10k , 100 * ( 1.0 * n10k / Scaffold_Number ) );
	fprintf ( fo, "scaffolds>100K\t%ld\t%.2f%\n", n100k, 100 * ( 1.0 * n100k / Scaffold_Number ) );
	fprintf ( fo, "scaffolds>1M  \t%ld\t%.2f%\n", n1m   , 100 * ( 1.0 * n1m / Scaffold_Number ) );
	fprintf ( fo, "\n" );
	fprintf ( fo, "Nucleotide_A\t%lld\t%.2f%\n", A_num_all, 100 * ( 1.0 * A_num_all / Size_includeN ) );
	fprintf ( fo, "Nucleotide_C\t%lld\t%.2f%\n", C_num_all, 100 * ( 1.0 * C_num_all / Size_includeN ) );
	fprintf ( fo, "Nucleotide_G\t%lld\t%.2f%\n", G_num_all, 100 * ( 1.0 * G_num_all / Size_includeN ) );
	fprintf ( fo, "Nucleotide_T\t%lld\t%.2f%\n", T_num_all, 100 * ( 1.0 * T_num_all / Size_includeN ) );
	fprintf ( fo, "GapContent_N\t%lld\t%.2f%\n", N_num_all, 100 * ( 1.0 * N_num_all / Size_includeN ) );
	fprintf ( fo, "Non_ACGTN\t%lld\t%.2f%\n", Non_ACGTN_all, 100 * ( 1.0 * Non_ACGTN_all / Size_includeN ) );
	fprintf ( fo, "GC_Content\t%.2f%\t\t(G+C)/(A+C+G+T)\n", 100 * ( 1.0 * ( G_num_all + C_num_all ) / ( A_num_all + C_num_all + G_num_all + T_num_all ) ) );
	fprintf ( fo, "\n" );

	for ( k = 0; k < 10; k++ )
		{ flag[k] = 0; }

	for ( k = Scaffold_Number - 1; k >= 0; k-- )
	{
		Sum = Sum + Size_Seq[k];

		if ( ( Sum >= Size_includeN * 0.1 ) && ( Sum < Size_includeN * 0.2 ) && ( flag[1] == 0 ) )
		{
			fprintf ( fo, "N10\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[1] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.2 ) && ( Sum < Size_includeN * 0.3 ) && ( flag[2] == 0 ) )
		{
			fprintf ( fo, "N20\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[2] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.3 ) && ( Sum < Size_includeN * 0.4 ) && ( flag[3] == 0 ) )
		{
			fprintf ( fo, "N30\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[3] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.4 ) && ( Sum < Size_includeN * 0.5 ) && ( flag[4] == 0 ) )
		{
			fprintf ( fo, "N40\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[4] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.5 ) && ( Sum < Size_includeN * 0.6 ) && ( flag[5] == 0 ) )
		{
			fprintf ( fo, "N50\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[5] = 1;
			N50 = Size_Seq[k];
		}
		else if ( ( Sum >= Size_includeN * 0.6 ) && ( Sum < Size_includeN * 0.7 ) && ( flag[6] == 0 ) )
		{
			fprintf ( fo, "N60\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[6] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.7 ) && ( Sum < Size_includeN * 0.8 ) && ( flag[7] == 0 ) )
		{
			fprintf ( fo, "N70\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[7] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.8 ) && ( Sum < Size_includeN * 0.9 ) && ( flag[8] == 0 ) )
		{
			fprintf ( fo, "N80\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[8] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.9 ) && ( flag[9] == 0 ) )
		{
			fprintf ( fo, "N90\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[9] = 1;
		}

		if ( ( Sum >= known_genome_size * 0.5 )  && ( flag_known == 0 ) )
		{
			N50_known = Size_Seq[k];
			Num_N50_known = Scaffold_Number - k;
			flag_known = 1;
		}
	}

	if ( flag[5] == 0 )
	{
		Sum = 0;

		for ( k = Scaffold_Number - 1; k >= 0; k-- )
		{
			Sum = Sum + Size_Seq[k];

			if ( Sum >= Size_includeN * 0.5 )
			{
				fprintf ( fo, "N50\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
				break;
			}
		}
	}

	fprintf ( fo, "\n" );

	if ( known_genome_size )
	{
		fprintf ( fo, "NG50\t%ld\t%ld\n", N50_known, Num_N50_known );
		fprintf ( fo, "N50_scaffold-NG50_scaffold_length_difference\t%ld\n", abs ( N50 - N50_known ) );
	}
	else
	{
		fprintf ( fo, "NG50\tNaN\tNaN\n" );
		fprintf ( fo, "N50_scaffold-NG50_scaffold_length_difference\tNaN\n" );
	}

	fprintf ( fo, "\n" );
	free ( A_num );
	free ( C_num );
	free ( G_num );
	free ( T_num );
	free ( N_num );
	free ( Non_ACGTN );
	free ( Singleton_Seq );
	free ( Size_Seq );
	Scaffold_Number_Scaf = Scaffold_Number;
	Singleton_Number_Scaf = Singleton_Number;
	/*********************** Contig ******************************/
	fprintf ( fo, "<-- Information for assembly Contig '%s.contig'.(cut_off_length < 100bp) -->\n\n", graphfile );
	cut_off_len = 0;
	Scaffold_Number = 0;
	Singleton_Number = 0;
	Singleton_Seq  = ( long * ) malloc ( sizeof ( long ) );
	A_num_all = 0;
	A_num  = ( long * ) malloc ( sizeof ( long ) );
	C_num_all = 0;
	C_num  = ( long * ) malloc ( sizeof ( long ) );
	G_num_all = 0;
	G_num  = ( long * ) malloc ( sizeof ( long ) );
	T_num_all = 0;
	T_num  = ( long * ) malloc ( sizeof ( long ) );
	N_num_all = 0;
	N_num  = ( long * ) malloc ( sizeof ( long ) );
	Non_ACGTN_all = 0;
	Non_ACGTN  = ( long * ) malloc ( sizeof ( long ) );
	Size_includeN = 0;
	Size_Seq        = ( long * ) malloc ( sizeof ( long ) );
	Sum = 0;
	n100 = 0;
	n500 = 0;
	n1k  = 0;
	n10k = 0;
	n100k = 0;
	n1m  = 0;
	N50 = 0;
	N50_known = 0;
	Num_N50_known = 0;
	flag_known = 0;
	cut_off_len = len_cut;
	A_num[Scaffold_Number] = 0;
	C_num[Scaffold_Number] = 0;
	G_num[Scaffold_Number] = 0;
	T_num[Scaffold_Number] = 0;
	N_num[Scaffold_Number] = 0;
	Non_ACGTN[Scaffold_Number] = 0;
	Size_Seq[Scaffold_Number] = 0;
	Singleton_Seq[Scaffold_Number] = 0;
	Nucleotide = fgetc ( fp2 );

	while ( Nucleotide != EOF )
	{
		if ( Nucleotide == '>' )
		{
			if ( ( Scaffold_Number > 0 ) && ( Size_Seq[Scaffold_Number - 1] < cut_off_len ) )
			{
				A_num_all = A_num_all - A_num[Scaffold_Number - 1];
				C_num_all = C_num_all - C_num[Scaffold_Number - 1];
				G_num_all = G_num_all - G_num[Scaffold_Number - 1];
				T_num_all = T_num_all - T_num[Scaffold_Number - 1];
				N_num_all = N_num_all - N_num[Scaffold_Number - 1];
				Non_ACGTN_all = Non_ACGTN_all - Non_ACGTN[Scaffold_Number - 1];
				Size_includeN = Size_includeN - Size_Seq[Scaffold_Number - 1];
				Singleton_Number = Singleton_Number - Singleton_Seq[Scaffold_Number - 1];
				Scaffold_Number = Scaffold_Number - 1;
			}
			else
			{
				Size_Seq = ( long * ) realloc ( Size_Seq, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				A_num    = ( long * ) realloc ( A_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				C_num    = ( long * ) realloc ( C_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				G_num    = ( long * ) realloc ( G_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				T_num    = ( long * ) realloc ( T_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				N_num    = ( long * ) realloc ( N_num, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				Non_ACGTN = ( long * ) realloc ( Non_ACGTN, ( Scaffold_Number + 2 ) * sizeof ( long ) );
				Singleton_Seq = ( long * ) realloc ( Singleton_Seq, ( Scaffold_Number + 2 ) * sizeof ( long ) );
			}

			Scaffold_Number++;
			A_num[Scaffold_Number - 1] = 0;
			C_num[Scaffold_Number - 1] = 0;
			G_num[Scaffold_Number - 1] = 0;
			T_num[Scaffold_Number - 1] = 0;
			N_num[Scaffold_Number - 1] = 0;
			Non_ACGTN[Scaffold_Number - 1] = 0;
			Size_Seq[Scaffold_Number - 1] = 0;
			Singleton_Seq[Scaffold_Number - 1] = 0;
			Nucleotide = fgetc ( fp2 );

			if ( Nucleotide == 'C' )
			{
				Singleton_Number++;
				Singleton_Seq[Scaffold_Number - 1] ++;
			}

			fgets ( buf, 4000, fp2 );
		}
		else if ( ( Nucleotide == 'N' ) || ( Nucleotide == 'n' ) )
		{
			N_num[Scaffold_Number - 1] ++;
			N_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'A' ) || ( Nucleotide == 'a' ) )
		{
			A_num[Scaffold_Number - 1] ++;
			A_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'C' ) || ( Nucleotide == 'c' ) )
		{
			C_num[Scaffold_Number - 1] ++;
			C_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'G' ) || ( Nucleotide == 'g' ) )
		{
			G_num[Scaffold_Number - 1] ++;
			G_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else if ( ( Nucleotide == 'T' ) || ( Nucleotide == 't' ) )
		{
			T_num[Scaffold_Number - 1] ++;
			T_num_all++;
			Size_Seq[Scaffold_Number - 1] ++;
			Size_includeN++;
		}
		else
		{
			if ( ( Nucleotide != '\n' ) && ( Nucleotide != '\r' ) )
			{
				Non_ACGTN[Scaffold_Number - 1] ++;
				Non_ACGTN_all++;
				Size_Seq[Scaffold_Number - 1] ++;
				Size_includeN++;
			}
		}

		Nucleotide = fgetc ( fp2 );
	}

	if ( Size_Seq[Scaffold_Number - 1] < cut_off_len )
	{
		A_num_all = A_num_all - A_num[Scaffold_Number - 1];
		C_num_all = C_num_all - C_num[Scaffold_Number - 1];
		G_num_all = G_num_all - G_num[Scaffold_Number - 1];
		T_num_all = T_num_all - T_num[Scaffold_Number - 1];
		N_num_all = N_num_all - N_num[Scaffold_Number - 1];
		Non_ACGTN_all = Non_ACGTN_all - Non_ACGTN[Scaffold_Number - 1];
		Size_includeN = Size_includeN - Size_Seq[Scaffold_Number - 1];
		Singleton_Number = Singleton_Number - Singleton_Seq[Scaffold_Number - 1];
		Scaffold_Number = Scaffold_Number - 1;
	}

	qsort ( Size_Seq, Scaffold_Number, sizeof ( Size_Seq[0] ), cmp_int );
	fprintf ( fo, "Size_includeN\t%lld\n", Size_includeN );
	fprintf ( fo, "Size_withoutN\t%lld\n", Size_includeN - N_num_all );
	fprintf ( fo, "Contig_Num\t%ld\n", Scaffold_Number );
	fprintf ( fo, "Mean_Size\t%lld\n", Size_includeN / Scaffold_Number );
	fprintf ( fo, "Median_Size\t%ld\n", Size_Seq[ ( Scaffold_Number + 1 ) / 2 - 1] );
	fprintf ( fo, "Longest_Seq\t%ld\n", Size_Seq[Scaffold_Number - 1] );
	fprintf ( fo, "Shortest_Seq\t%ld\n", Size_Seq[0] );
	fprintf ( fo, "\n" );

	for ( k = 0; k < Scaffold_Number; k++ )
	{
		if ( Size_Seq[k] > 100 )
		{
			n100++;
		}

		if ( Size_Seq[k] > 500 )
		{
			n500++;
		}

		if ( Size_Seq[k] > 1000 )
		{
			n1k++;
		}

		if ( Size_Seq[k] > 10000 )
		{
			n10k++;
		}

		if ( Size_Seq[k] > 100000 )
		{
			n100k++;
		}

		if ( Size_Seq[k] > 1000000 )
		{
			n1m++;
		}
	}

	fprintf ( fo, "Contig>100 \t%ld\t%.2f%\n", n100 , 100 * ( 1.0 * n100 / Scaffold_Number ) );
	fprintf ( fo, "Contig>500 \t%ld\t%.2f%\n", n500 , 100 * ( 1.0 * n500 / Scaffold_Number ) );
	fprintf ( fo, "Contig>1K  \t%ld\t%.2f%\n", n1k      , 100 * ( 1.0 * n1k / Scaffold_Number ) );
	fprintf ( fo, "Contig>10K \t%ld\t%.2f%\n", n10k , 100 * ( 1.0 * n10k / Scaffold_Number ) );
	fprintf ( fo, "Contig>100K\t%ld\t%.2f%\n", n100k, 100 * ( 1.0 * n100k / Scaffold_Number ) );
	fprintf ( fo, "Contig>1M  \t%ld\t%.2f%\n", n1m      , 100 * ( 1.0 * n1m / Scaffold_Number ) );
	fprintf ( fo, "\n" );
	fprintf ( fo, "Nucleotide_A\t%lld\t%.2f%\n", A_num_all, 100 * ( 1.0 * A_num_all / Size_includeN ) );
	fprintf ( fo, "Nucleotide_C\t%lld\t%.2f%\n", C_num_all, 100 * ( 1.0 * C_num_all / Size_includeN ) );
	fprintf ( fo, "Nucleotide_G\t%lld\t%.2f%\n", G_num_all, 100 * ( 1.0 * G_num_all / Size_includeN ) );
	fprintf ( fo, "Nucleotide_T\t%lld\t%.2f%\n", T_num_all, 100 * ( 1.0 * T_num_all / Size_includeN ) );
	fprintf ( fo, "GapContent_N\t%lld\t%.2f%\n", N_num_all, 100 * ( 1.0 * N_num_all / Size_includeN ) );
	fprintf ( fo, "Non_ACGTN\t%lld\t%.2f%\n", Non_ACGTN_all, 100 * ( 1.0 * Non_ACGTN_all / Size_includeN ) );
	fprintf ( fo, "GC_Content\t%.2f%\t\t(G+C)/(A+C+G+T)\n", 100 * ( 1.0 * ( G_num_all + C_num_all ) / ( A_num_all + C_num_all + G_num_all + T_num_all ) ) );
	fprintf ( fo, "\n" );

	for ( k = 0; k < 10; k++ )
		{ flag[k] = 0; }

	for ( k = Scaffold_Number - 1; k >= 0; k-- )
	{
		Sum = Sum + Size_Seq[k];

		if ( ( Sum >= Size_includeN * 0.1 ) && ( Sum < Size_includeN * 0.2 ) && ( flag[1] == 0 ) )
		{
			fprintf ( fo, "N10\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[1] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.2 ) && ( Sum < Size_includeN * 0.3 ) && ( flag[2] == 0 ) )
		{
			fprintf ( fo, "N20\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[2] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.3 ) && ( Sum < Size_includeN * 0.4 ) && ( flag[3] == 0 ) )
		{
			fprintf ( fo, "N30\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[3] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.4 ) && ( Sum < Size_includeN * 0.5 ) && ( flag[4] == 0 ) )
		{
			fprintf ( fo, "N40\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[4] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.5 ) && ( Sum < Size_includeN * 0.6 ) && ( flag[5] == 0 ) )
		{
			fprintf ( fo, "N50\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[5] = 1;
			N50 = Size_Seq[k];
		}
		else if ( ( Sum >= Size_includeN * 0.6 ) && ( Sum < Size_includeN * 0.7 ) && ( flag[6] == 0 ) )
		{
			fprintf ( fo, "N60\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[6] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.7 ) && ( Sum < Size_includeN * 0.8 ) && ( flag[7] == 0 ) )
		{
			fprintf ( fo, "N70\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[7] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.8 ) && ( Sum < Size_includeN * 0.9 ) && ( flag[8] == 0 ) )
		{
			fprintf ( fo, "N80\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[8] = 1;
		}
		else if ( ( Sum >= Size_includeN * 0.9 ) && ( flag[9] == 0 ) )
		{
			fprintf ( fo, "N90\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
			flag[9] = 1;
		}

		if ( ( Sum >= known_genome_size * 0.5 )  && ( flag_known == 0 ) )
		{
			N50_known = Size_Seq[k];
			Num_N50_known = Scaffold_Number - k;
			flag_known = 1;
		}
	}

	if ( flag[5] == 0 )
	{
		Sum = 0;

		for ( k = Scaffold_Number - 1; k >= 0; k-- )
		{
			Sum = Sum + Size_Seq[k];

			if ( Sum >= Size_includeN * 0.5 )
			{
				fprintf ( fo, "N50\t%ld\t%ld\n", Size_Seq[k], Scaffold_Number - k );
				break;
			}
		}
	}

	fprintf ( fo, "\n" );

	if ( known_genome_size )
	{
		fprintf ( fo, "NG50\t%ld\t%ld\n", N50_known, Num_N50_known );
		fprintf ( fo, "N50_contig-NG50_contig_length_difference\t%ld\n", abs ( N50 - N50_known ) );
	}
	else
	{
		fprintf ( fo, "NG50\tNaN\tNaN\n" );
		fprintf ( fo, "N50_contig-NG50_contig_length_difference\tNaN\n" );
	}

	fprintf ( fo, "\n" );
	free ( A_num );
	free ( C_num );
	free ( G_num );
	free ( T_num );
	free ( N_num );
	free ( Non_ACGTN );
	free ( Singleton_Seq );
	free ( Size_Seq );
	Scaffold_Number_Contig = Scaffold_Number;
	fprintf ( fo, "Number_of_contigs_in_scaffolds\t%ld\n", Scaffold_Number_Contig - Singleton_Number_Scaf );
	fprintf ( fo, "Number_of_contigs_not_in_scaffolds(Singleton)\t%ld\n", Singleton_Number_Scaf );
	fprintf ( fo, "Average_number_of_contigs_per_scaffold\t%.1f\n", 1.0 * ( Scaffold_Number_Contig - Singleton_Number_Scaf ) / ( Scaffold_Number_Scaf - Singleton_Number_Scaf ) );
	fprintf ( fo, "\n" );
	fclose ( fp );
	fclose ( fp2 );
	fclose ( fo );
}

/*************************************************
Function:
    getCnt
Description:
    Gets connection between two contigs.
Input:
    1. from_c:      left contig
    2. to_c:             right contig
Output:
    None.
Return:
    Connection between two contigs.
*************************************************/
CONNECT * getCnt ( unsigned int from_c, unsigned int to_c )
{
	CONNECT * pcnt;
	pcnt = contig_array[from_c].downwardConnect;

	while ( pcnt )
	{
		if ( pcnt->contigID == to_c )
			{ return pcnt; }

		pcnt = pcnt->next;
	}

	return pcnt;
}

/*************************************************
Function:
    allConnect
Description:
    Counts contig's downstrem connection number if this contig does
    NOT have neigher upstream nor downstream connection in scaffold.
Input:
    1. ctg:         contig
    2. preCNT:       contig's upstream connection
Output:
    None.
Return:
    1 if contig had both upstream and downstream connections in
       scaffold.
    Contig's downstream connection number otherwise.
*************************************************/
int allConnect ( unsigned int ctg, CONNECT * preCNT )
{
	if ( preCNT && preCNT->nextInScaf )
		{ return 1; }

	CONNECT * cn_temp;
	int count = 0;

	if ( !contig_array[ctg].downwardConnect )
		{ return count; }

	cn_temp = contig_array[ctg].downwardConnect;

	while ( cn_temp )
	{
		count++;
		cn_temp = cn_temp->next;
	}

	return count;
}

/*************************************************
Function:
    get_ctg_score2
Description:
    Calculates a score to indicate the strengh of connections
    between a contig and other contigs in scaffold.
Input:
    1. pos:             contig index
    2. tempCounter:     maximum contig number in scaffold
Output:
    None.
Return:
    Calculated score.
*************************************************/
int get_ctg_score2 ( int pos, int tempCounter )
{
	int id, innum, outnum, in = 0, out = 0, i, currId;
	CONNECT * dh_cnt;
	id = * ( unsigned int * ) darrayGet ( tempArray, pos );
	outnum = allConnect ( id, NULL );
	innum = allConnect ( getTwinCtg ( id ), NULL );
	int outlen = 0, inlen = 0, outcut, incut;

	if ( contig_array[id].downwardConnect )
		{ outlen = contig_array[id].downwardConnect->gapLen; }

	if ( contig_array[getTwinCtg ( id )].downwardConnect )
		{ inlen = contig_array[getTwinCtg ( id )].downwardConnect->gapLen; }

	outcut = outlen / overlaplen;
	incut = inlen / overlaplen;

	if ( outcut > outnum * 10 || outcut < outnum * 2 )
		{ outcut = outnum * 10; }

	if ( incut > innum * 10 || incut < innum * 2 )
		{ incut = innum * 10; }

	int start = pos - incut > 0 ? pos - incut : 0;
	int end = pos + outcut < tempCounter ? pos + outcut : tempCounter;

	for ( i = start; i < end; i++ )
	{
		if ( i == pos )
			{ continue; }

		currId = * ( unsigned int * ) darrayGet ( tempArray, i );

		if ( i < pos )
		{
			dh_cnt = getCnt ( currId, id );

			if ( dh_cnt && dh_cnt->weight > 0 )
				{ in++; }
		}

		if ( i > pos )
		{
			dh_cnt = getCnt ( id, currId );

			if ( dh_cnt && dh_cnt->weight > 0 )
				{ out++; }
		}
	}

	if ( innum > pos )
		{ innum = pos; }

	if ( outnum > tempCounter - pos - 1 )
		{ outnum = tempCounter - pos - 1; }

	if ( innum > 0 && outnum > 0 )
	{
		if ( pos < tempCounter - 1 && pos > 0 )
		{
			if ( outnum < 3 && out == 0 )
			{
				if ( in > 0 )
					{ return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 ); }
			}

			if ( innum < 3 && in == 0 )
			{
				if ( out > 0 )
					{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }
			}

			if ( in == 0 || out == 0 )
				{ return 0; }

			return ( int ) ( ( ( double ) ( in * out ) / ( double ) ( innum * outnum ) ) * 100 );
		}

		if ( pos == 0 )
			{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }

		if ( pos == tempCounter - 1 )
			{ return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 ); }
	}
	else if ( innum > 0 )
	{
		if ( pos == tempCounter - 1 && in == 1 && innum > 5 )
			{ return 0; }

		return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 );
	}
	else if ( outnum > 0 )
	{
		if ( pos == 0 && out == 1 && outnum > 5 )
			{ return 0; }

		return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 );
	}

	return 0;
}

/*************************************************
Function:
    get_ctg_score2
Description:
    Calculates a score to indicate the strengh of connections
    between a contig and other contigs in scaffold.
Input:
    1. pos:             contig index
    2. tempCounter:     maximum contig number in scaffold
Output:
    None.
Return:
    Calculated score.
*************************************************/
int get_ctg_score ( int pos, int num3, int num5, int flag )
{
	int i = 0, in = 0, out = 0, innum = 0, outnum = 0, threeid;
	CONNECT * dh_cnt;
	int id;

	if ( flag == 0 )
	{
		id = * ( unsigned int * ) darrayGet ( scaf3, pos );
		int end = pos > 100 ? pos - 100 : 0;
		int start = pos + 100 < num3 ? pos + 100 : num3;
		in = 0, out = 0;

		for ( i = start; i >= end; i-- )
		{
			threeid = * ( unsigned int * ) darrayGet ( scaf3, i );

			if ( threeid == id )
			{
				pos = i;
				continue;
			}

			dh_cnt = getCnt ( id, threeid );

			if ( dh_cnt && dh_cnt->weight > 0 )
			{
				out++;
			}

			dh_cnt = getCnt ( threeid, id );

			if ( dh_cnt && dh_cnt->weight > 0 )
				{ in++; }
		}

		outnum = allConnect ( id, NULL );
		innum = allConnect ( getTwinCtg ( id ), NULL );
		int num5_check = 0;

		if ( pos - end < outnum )
		{
			for ( i = 0; i < num5; i++ )
			{
				num5_check++;
				threeid = * ( unsigned int * ) darrayGet ( scaf5, i );
				dh_cnt = getCnt ( id, threeid );

				if ( dh_cnt && dh_cnt->weight > 0 )
					{ out++; }

				if ( num5_check == outnum )
					{ break; }
			}
		}

		if ( pos - end + num5_check < outnum )
			{ outnum = pos - end + num5_check; }

		if ( start - pos < innum )
			{ innum = start - pos; }

		if ( innum > 0 && outnum > 0 )
		{
			if ( pos != num3 && pos > 0 )
			{
				if ( outnum < 5 && out == 0 )
				{
					if ( innum > 0 )
						{ return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 ); }
				}

				if ( innum < 5 && in == 0 )
				{
					if ( outnum > 0 )
						{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }
				}

				if ( in == 0 || out == 0 )
					{ return 0; }

				return ( int ) ( ( ( double ) ( in * out ) / ( double ) ( innum * outnum ) ) * 100 );
			}

			if ( pos == num3 )
				{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }

			if ( pos == 0 )
				{ return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 ); }
		}
		else if ( innum > 0 )
		{
			return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 );
		}
		else if ( outnum > 0 )
		{
			if ( pos == num3 && out == 1 && outnum > 5 )
				{ return 0; }

			return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 );
		}

		return 0;
	}
	else
	{
		id = * ( unsigned int * ) darrayGet ( scaf5, pos );
		int start = pos > 100 ? pos - 100 : 0;
		int end = pos + 100 < num5 ? pos + 100 : num5;
		in = 0, out = 0;

		for ( i = start; i < end; i++ )
		{
			threeid = * ( unsigned int * ) darrayGet ( scaf5, i );

			if ( threeid == id )
			{
				pos = i;
				continue;
			}

			dh_cnt = getCnt ( id, threeid );

			if ( dh_cnt && dh_cnt->weight > 0 )
				{ out++; }

			dh_cnt = getCnt ( threeid, id );

			if ( dh_cnt && dh_cnt->weight > 0 )
				{ in++; }
		}

		outnum = allConnect ( id, NULL );
		innum = allConnect ( getTwinCtg ( id ), NULL );
		int num3_check = 0;

		if ( pos - start < innum )
		{
			for ( i = 0; i < num3; i++ )
			{
				num3_check++;
				threeid = * ( unsigned int * ) darrayGet ( scaf3, i );
				dh_cnt = getCnt ( threeid, id );

				if ( dh_cnt && dh_cnt->weight > 0 )
					{ in++; }

				if ( num3_check == innum )
					{ break; }
			}
		}

		if ( pos - start + num3_check < innum )
			{ innum = pos - start + num3_check; }

		if ( end - pos - 1 < outnum )
			{ outnum = end - pos - 1; }

		if ( innum > 0 && outnum > 0 )
		{
			if ( pos != num5 - 1 && pos > 0 )
			{
				if ( outnum < 5 && out == 0 && innum > 0 )
					{ return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 ); }

				if ( innum < 5 && in == 0 && outnum > 0 )
					{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }

				if ( in == 0 || out == 0 )
					{ return 0; }

				return ( int ) ( ( ( double ) ( in * out ) / ( double ) ( innum * outnum ) ) * 100 );
			}

			if ( pos == 0 )
				{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }

			if ( pos == num5 - 1 )
				{ return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 ); }
		}
		else if ( innum > 0 )
		{
			if ( pos == num5 - 1 && in == 1 && innum > 5 )
				{ return 0; }

			return ( int ) ( ( ( double ) ( in ) / ( double ) ( innum ) ) * 100 );
		}
		else if ( outnum > 0 )
			{ return ( int ) ( ( ( double ) ( out ) / ( double ) ( outnum ) ) * 100 ); }

		return 0;
	}

	return 0;
}

/*************************************************
Function:
    scaffolding
Description:
    Masks contigs having low connection score in scaffold. Fills gaps
    by using ARC information. Outputs scaffold structure and filled
    gaps' information.
Input:
    1. len_cut:     length cutoff
    2. outfile:         prefix of graph
Output:
    None.
Return:
    None.
*************************************************/
void scaffolding ( unsigned int len_cut, char * outfile )
{
	unsigned int prev_ctg, ctg, bal_ctg, *length_array, count = 0, num_lctg = 0, *score_array;
	unsigned int i, max_steps = 5;
	int num5, num3, j, len, flag, num_route, gap_c = 0;
	int tempCounter;
	short gap = 0;
	long long sum = 0, N50, N90;
	FILE * fp, *fo = NULL;
	char name[256];
	CONNECT * cnt, *prevCNT, *nextCnt, *dh_cnt;
	boolean excep, weak;
	weakCounter = 0;
	so_far = ( unsigned int * ) ckalloc ( max_n_routes * sizeof ( unsigned int ) );
	found_routes = ( unsigned int ** ) ckalloc ( max_n_routes * sizeof ( unsigned int * ) );

	for ( j = 0; j < max_n_routes; j++ )
		{ found_routes[j] = ( unsigned int * ) ckalloc ( max_steps * sizeof ( unsigned int ) ); }

	length_array = ( unsigned int * ) ckalloc ( ( num_ctg + 1 ) * sizeof ( unsigned int ) );

	//use length_array to change info in index_array
	for ( i = 1; i <= num_ctg; i++ )
		{ length_array[i] = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( index_array[i] > 0 )
			{ length_array[index_array[i]] = i; }
	}

	for ( i = 1; i <= num_ctg; i++ )
		{ index_array[i] = length_array[i]; }

	orig2new = 0;
	sprintf ( name, "%s.scaf", outfile );
	fp = ckopen ( name, "w" );
	sprintf ( name, "%s.scaf_gap", outfile );
	fo = ckopen ( name, "w" );
	scaf3 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	scaf5 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	gap3 = ( DARRAY * ) createDarray ( 1000, sizeof ( int ) );
	gap5 = ( DARRAY * ) createDarray ( 1000, sizeof ( int ) );
	tempArray = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );

	for ( i = 1; i <= num_ctg; i++ )
		{ contig_array[i].flag = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].length + ( unsigned int ) overlaplen >= len_cut )
			{ num_lctg++; }
		else
			{ continue; }

		if ( contig_array[i].flag || contig_array[i].mask || !contig_array[i].downwardConnect || !validConnect ( i, NULL ) )
			{ continue; }

		num5 = num3 = 0;
		ctg = i;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = i;
		contig_array[i].flag = 1;
		bal_ctg = getTwinCtg ( ctg );
		contig_array[bal_ctg].flag = 1;
		len = contig_array[i].length;
		prevCNT = NULL;
		cnt = getNextContig ( ctg, prevCNT, &excep );

		while ( cnt )
		{
			nextCnt = getNextContig ( cnt->contigID, cnt, &excep );

			if ( excep && prevCNT )
				{ fprintf ( stderr, "scaffolding: exception --- prev cnt from %u\n", prevCNT->contigID ); }

			if ( nextCnt && nextCnt->used )
				{ break; }

			setConnectUsed ( ctg, cnt->contigID, 1 );
			* ( int * ) darrayPut ( gap5, num5 - 1 ) = cnt->gapLen;
			ctg = cnt->contigID;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg;
			len += cnt->gapLen + contig_array[ctg].length;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		ctg = getTwinCtg ( i );

		if ( num5 >= 2 )
			{ prevCNT = checkConnect ( getTwinCtg ( * ( unsigned int * ) darrayGet ( scaf5, 1 ) ), ctg ); }
		else
			{ prevCNT = NULL; }

		cnt = getNextContig ( ctg, prevCNT, &excep );

		while ( cnt )
		{
			nextCnt = getNextContig ( cnt->contigID, cnt, &excep );

			if ( excep && prevCNT )
				{ fprintf ( stderr, "scaffolding: exception -- prev cnt from %u\n", prevCNT->contigID ); }

			if ( nextCnt && nextCnt->used )
				{ break; }

			setConnectUsed ( ctg, cnt->contigID, 1 );
			ctg = cnt->contigID;
			len += cnt->gapLen + contig_array[ctg].length;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			* ( int * ) darrayPut ( gap3, num3 ) = cnt->gapLen;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = bal_ctg;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		if ( num5 + num3 == 1 )
		{
			contig_array[i].flag = 0;
			continue;
		}

		len += overlaplen;
		sum += len;
		length_array[count++] = len;

		if ( num5 + num3 < 1 )
		{
			fprintf ( stderr, "no scaffold created for contig %d\n", i );
			continue;
		}

		tempCounter = 0;

		for ( j = num3 - 1; j >= 0; j-- )
			{ * ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) = * ( unsigned int * ) darrayGet ( scaf3, j ); }

		for ( j = 0; j < num5; j++ )
			{ * ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) = * ( unsigned int * ) darrayGet ( scaf5, j ); }

		score_array = ( unsigned int * ) ckalloc ( tempCounter * sizeof ( unsigned int ) );
		int now_cnt_weight, curr_ctg_score, pre_score = -1, prev_id = 0, mask_num = 0, score_count = 0, prev_p = 0;

		for ( j = 0; j < tempCounter; j++ )
		{
			int currId = * ( unsigned int * ) darrayGet ( tempArray, j );
			dh_cnt = getCntBetween ( prev_id, currId );

			if ( dh_cnt )
				{ now_cnt_weight = dh_cnt ->weight; }
			else
				{ now_cnt_weight = 0; }

			curr_ctg_score = get_ctg_score2 ( j, tempCounter );

			if ( prev_id == 0 )
			{
				pre_score = curr_ctg_score;
				prev_id = currId;
				prev_p = j;
				continue;
			}

			if ( score_mask )
			{
				if ( now_cnt_weight == 0 && j > 0 && j < tempCounter - 1 )
				{
					if ( pre_score == 0 )
					{
						* ( unsigned int * ) darrayPut ( tempArray, prev_p ) = 0;
						contig_array[prev_id].flag = 0;
						contig_array[getTwinCtg ( prev_id )].flag = 0;
						mask_num++;
					}
					else if ( curr_ctg_score == 0 )
					{
						* ( unsigned int * ) darrayPut ( tempArray, j ) = 0;
						contig_array[currId].flag = 0;
						contig_array[getTwinCtg ( currId )].flag = 0;
						mask_num++;

						if ( j < tempCounter - 1 )
							{ continue; }
					}
				}

				if ( abs ( prev_id - currId ) <= 2 && now_cnt_weight == 0
				        && * ( unsigned int * ) darrayGet ( tempArray, prev_p ) != 0 && * ( unsigned int * ) darrayGet ( tempArray, j ) != 0 )
				{
					mask_num++;

					if ( contig_array[prev_id].cvg < contig_array[currId].cvg )
					{
						* ( unsigned int * ) darrayPut ( tempArray, prev_p ) = 0;
						contig_array[prev_id].flag = 0;
						contig_array[getTwinCtg ( prev_id )].flag = 0;
					}
					else
					{
						* ( unsigned int * ) darrayPut ( tempArray, j ) = 0;
						contig_array[currId].flag = 0;
						contig_array[getTwinCtg ( currId )].flag = 0;

						if ( j < tempCounter - 1 )
							{ continue; }
					}
				}
			}

			if ( * ( unsigned int * ) darrayGet ( tempArray, prev_p ) != 0 )
				{ score_array[score_count++] = pre_score; }

			pre_score = curr_ctg_score;
			prev_id = currId;
			prev_p = j;
		}

		if ( * ( unsigned int * ) darrayGet ( tempArray, prev_p ) != 0 )
			{ score_array[score_count++] = pre_score; }

		if ( score_mask == 1 && ( ( num3 + num5 > 5 && score_count < 2 ) || score_count == 1 ) )
		{
			free ( ( void * ) score_array );
			--count;
			sum -= len;

			for ( j = 0; j < num3; j++ )
			{
				ctg = * ( unsigned int * ) darrayGet ( scaf3, j );
				contig_array[ctg].flag = 0;
				contig_array[getTwinCtg ( ctg )].flag = 0;
			}

			for ( j = 0; j < num5; j++ )
			{
				ctg = * ( unsigned int * ) darrayGet ( scaf5, j );
				contig_array[ctg].flag = 0;
				contig_array[getTwinCtg ( ctg )].flag = 0;
			}

			continue;
		}

		fprintf ( fp, ">scaffold%d %d %d %d\n", count, score_count, len, num3 + num5, mask_num );
		fprintf ( fo, ">scaffold%d %d %d %d\n", count, score_count, len, num3 + num5, mask_num );
		len = prev_ctg = 0;
		tempCounter = 0, score_count = 0;

		for ( j = num3 - 1; j >= 0; j-- )
		{
			int now_cnt_weigth = 0;
			int nextid, start = 0;
			int currId = * ( unsigned int * ) darrayGet ( scaf3, j );
			int tmpid = * ( unsigned int * ) darrayGet ( tempArray, tempCounter++ );

			if ( tmpid == 0 )
			{
				if ( j == num3 - 1 )
				{
					len = 0;
					continue;
				}

				len += contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + * ( int * ) darrayGet ( gap3, j );
				int tmpgap = contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + * ( int * ) darrayGet ( gap3, j );
				gap += tmpgap > 0 ? tmpgap : 0;
				continue;
			}

			if ( j > 0 )
			{
				nextid = * ( unsigned int * ) darrayGet ( tempArray, tempCounter + start );

				while ( nextid == 0 && tempCounter + start + 1 < num3 + num5 )
				{
					start++;
					nextid = * ( unsigned int * ) darrayGet ( tempArray, tempCounter + start );
				}
			}
			else
			{
				nextid = i;
			}

			CONNECT * dh_cnt = getCntBetween ( currId, nextid );

			if ( dh_cnt )
				{ now_cnt_weigth = dh_cnt->weight; }
			else
				{ now_cnt_weigth = 0; }

			curr_ctg_score = score_array[score_count++];

			if ( score_mask == 1 && curr_ctg_score == 0 && ( num3 + num5 > 2 )
			        && ( ( j == num3 - 1 && contig_array[nextid].length < 200 && num3 + num5 > 5 ) || ( now_cnt_weigth == 0 && j > 0 ) ) )
			{
				if ( j == num3 - 1 )
				{
					len = 0;
					continue;
				}

				len += contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + * ( int * ) darrayGet ( gap3, j );
				int tmpgap = contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + * ( int * ) darrayGet ( gap3, j );
				gap += tmpgap > 0 ? tmpgap : 0;
				continue;
			}

			if ( !isLargerThanTwin ( * ( unsigned int * ) darrayGet ( scaf3, j ) ) )
			{
				fprintf ( fp, "%-10d %-10d +   %d  %d  %d"
				          , index_array[* ( unsigned int * ) darrayGet ( scaf3, j )], len,
				          contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + overlaplen,
				          now_cnt_weigth, curr_ctg_score );
				weak = printCnts ( fp, * ( unsigned int * ) darrayGet ( scaf3, j ) );
			}
			else
			{
				fprintf ( fp, "%-10d %-10d -   %d  %d  %d"
				          , index_array[getTwinCtg ( * ( unsigned int * ) darrayGet ( scaf3, j ) )], len
				          , contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + overlaplen,
				          now_cnt_weigth, curr_ctg_score );
				weak = printCnts ( fp, * ( unsigned int * ) darrayGet ( scaf3, j ) );
			}

			if ( prev_ctg )
			{
				num_route = num_trace = 0;
				traceAlongArc ( * ( unsigned int * ) darrayGet ( scaf3, j ), prev_ctg, max_steps
				                , gap - ins_size_var, gap + ins_size_var, 0, 0, &num_route );

				if ( num_route == 1 )
				{
					output1gap ( fo, max_steps );
					gap_c++;
				}
			}

			fprintf ( fo, "%-10d %-10d\n", * ( unsigned int * ) darrayGet ( scaf3, j ), len );
			len += contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + * ( int * ) darrayGet ( gap3, j );
			prev_ctg = * ( unsigned int * ) darrayGet ( scaf3, j );
			gap = * ( int * ) darrayGet ( gap3, j ) > 0 ? * ( int * ) darrayGet ( gap3, j ) : 0;
		}

		for ( j = 0; j < num5; j++ )
		{
			int now_cnt_weigth = 0;
			int currId, nextid, start = 0;
			currId = * ( unsigned int * ) darrayGet ( scaf5, j );
			int tmpid = * ( unsigned int * ) darrayGet ( tempArray, tempCounter++ );

			if ( tmpid == 0 )
			{
				if ( j == num5 - 1 )
					{ continue; }

				len += contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length + * ( int * ) darrayGet ( gap5, j );
				int tmpgap = contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length + * ( int * ) darrayGet ( gap5, j );
				gap += tmpgap > 0 ? tmpgap : 0;
				continue;
			}

			if ( j < num5 - 1 )
			{
				nextid = * ( unsigned int * ) darrayGet ( tempArray, tempCounter + start );

				while ( nextid == 0 && tempCounter + start + 1 < num3 + num5 )
				{
					start ++;
					nextid = * ( unsigned int * ) darrayGet ( tempArray, tempCounter + start );
				}

				CONNECT * dh_cnt = getCntBetween ( currId, nextid );

				if ( dh_cnt )
					{ now_cnt_weigth = dh_cnt->weight; }
				else
					{ now_cnt_weigth = 0; }
			}

			curr_ctg_score = score_array [score_count++];

			if ( score_mask == 1 && curr_ctg_score == 0 && ( num3 + num5 > 2 )
			        && ( ( j == num5 - 1 && contig_array[* ( unsigned int * ) darrayGet ( scaf5, j - 1 )].length < 200 && num3 + num5 > 5 )
			             || ( j != num5 - 1 && now_cnt_weigth == 0 ) ) )
			{
				if ( j == num5 - 1 )
				{
					continue;
				}

				len += contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length + * ( int * ) darrayGet ( gap5, j );
				int tmpgap = contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length + * ( int * ) darrayGet ( gap5, j );
				gap += tmpgap > 0 ? tmpgap : 0;
				continue;
			}

			if ( !isLargerThanTwin ( * ( unsigned int * ) darrayGet ( scaf5, j ) ) )
			{
				fprintf ( fp, "%-10d %-10d +   %d  %d  %d"
				          , index_array[* ( unsigned int * ) darrayGet ( scaf5, j )], len
				          , contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length + overlaplen,
				          now_cnt_weigth, curr_ctg_score );
				weak = printCnts ( fp, * ( unsigned int * ) darrayGet ( scaf5, j ) );
			}
			else
			{
				fprintf ( fp, "%-10d %-10d -   %d  %d  %d"
				          , index_array[getTwinCtg ( * ( unsigned int * ) darrayGet ( scaf5, j ) )], len
				          , contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length + overlaplen,
				          now_cnt_weigth, curr_ctg_score );
				weak = printCnts ( fp, * ( unsigned int * ) darrayGet ( scaf5, j ) );
			}

			if ( prev_ctg )
			{
				num_route = num_trace = 0;
				traceAlongArc ( * ( unsigned int * ) darrayGet ( scaf5, j ), prev_ctg, max_steps
				                , gap - ins_size_var, gap + ins_size_var, 0, 0, &num_route );

				if ( num_route == 1 )
				{
					output1gap ( fo, max_steps );
					gap_c++;
				}
			}

			fprintf ( fo, "%-10d %-10d\n", * ( unsigned int * ) darrayGet ( scaf5, j ), len );

			if ( j < num5 - 1 )
			{
				len += contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length +
				       * ( int * ) darrayGet ( gap5, j );
				prev_ctg = * ( unsigned int * ) darrayGet ( scaf5, j );
				gap = * ( int * ) darrayGet ( gap5, j ) > 0 ? * ( int * ) darrayGet ( gap5, j ) : 0;
			}
		}

		free ( ( void * ) score_array );
	}

	freeDarray ( scaf3 );
	freeDarray ( scaf5 );
	freeDarray ( gap3 );
	freeDarray ( gap5 );
	freeDarray ( tempArray );
	fclose ( fp );
	fclose ( fo );
	fprintf ( stderr, "\nThe final rank\n" );

	if ( count == 0 )
	{
		fprintf ( stderr, "\n\nNo scaffold was constructed.\n\n" );
		free ( ( void * ) length_array );

		for ( j = 0; j < max_n_routes; j++ )
			{ free ( ( void * ) found_routes[j] ); }

		free ( ( void * ) found_routes );
		free ( ( void * ) so_far );
		return;
	}

	fprintf ( stderr, "\n*******************************\n" );
	fprintf ( stderr, " Scaffold number                  %d\n", count );
	fprintf ( stderr, " In-scaffold contig number        %u\n", num_lctg / 2 );
	fprintf ( stderr, " Total scaffold length            %lld\n", sum );
	fprintf ( stderr, " Average scaffold length          %lld\n", sum / count );
	fprintf ( stderr, " Filled gap number                %d\n", gap_c );

	//output singleton
	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].length + ( unsigned int ) overlaplen < len_cut || contig_array[i].flag )
			{ continue; }

		length_array[count++] = contig_array[i].length;
		sum += contig_array[i].length;

		if ( isSmallerThanTwin ( i ) )
			{ i++; }
	}

	long long total_len = sum;
	qsort ( length_array, count, sizeof ( length_array[0] ), cmp_int );
	N50 = sum * 0.5;
	N90 = sum * 0.9;
	int N50length = 0;
	int N90length = 0;
	sum = flag = 0;

	for ( j = count - 1; j >= 0; j-- )
	{
		sum += length_array[j];

		if ( !flag && sum >= N50 && N50length == 0 )
		{
			N50length = length_array[j];
			flag++;
		}

		if ( sum >= N90 && N90length == 0 )
		{
			N90length = length_array[j];
			break;
		}
	}

	fprintf ( stderr, " Longest scaffold                 %lld\n", length_array[count - 1] );
	fprintf ( stderr, " Scaffold and singleton number    %d\n", count );
	fprintf ( stderr, " Scaffold and singleton length    %lld\n", total_len );
	fprintf ( stderr, " Average length                   %d\n", total_len / count );
	fprintf ( stderr, " N50                              %d\n", N50length );
	fprintf ( stderr, " N90                              %d\n", N90length );
	fprintf ( stderr, " Weak points                      %d\n", weakCounter );
	fprintf ( stderr, "\n*******************************\n" );
	fflush ( stdout );
	free ( ( void * ) length_array );

	for ( j = 0; j < max_n_routes; j++ )
		{ free ( ( void * ) found_routes[j] ); }

	free ( ( void * ) found_routes );
	free ( ( void * ) so_far );
}

/*************************************************
Function:
    scaffold_count
Description:
    Makes statistic of scaffold till current rank.
Input:
    1. rank:            current rank number
    2. len_cut:     length cutoff
Output:
    None.
Return:
    None.
*************************************************/
void scaffold_count ( int rank, unsigned int len_cut )
{
	static DARRAY * scaf3, *scaf5;
	static DARRAY * gap3, *gap5;
	unsigned int prev_ctg, ctg, bal_ctg, *length_array, count = 0, num_lctg = 0;
	unsigned int i, max_steps = 5;
	int num5, num3, j, len, flag, num_route, gap_c = 0;
	short gap = 0;
	long long sum = 0, N50, N90;
	CONNECT * cnt, *prevCNT, *nextCnt;
	boolean excep;
	so_far = ( unsigned int * ) ckalloc ( max_n_routes * sizeof ( unsigned int ) );
	found_routes = ( unsigned int ** ) ckalloc ( max_n_routes * sizeof ( unsigned int * ) );

	for ( j = 0; j < max_n_routes; j++ )
		{ found_routes[j] = ( unsigned int * ) ckalloc ( max_steps * sizeof ( unsigned int ) ); }

	length_array = ( unsigned int * ) ckalloc ( ( num_ctg + 1 ) * sizeof ( unsigned int ) );

	//use length_array to change info in index_array
	for ( i = 1; i <= num_ctg; i++ )
		{ length_array[i] = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( index_array[i] > 0 )
			{ length_array[index_array[i]] = i; }
	}

	for ( i = 1; i <= num_ctg; i++ )
		{ index_array[i] = length_array[i]; }  //contig i with original index: index_array[i]

	orig2new = 0;
	scaf3 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	scaf5 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	gap3 = ( DARRAY * ) createDarray ( 1000, sizeof ( int ) );
	gap5 = ( DARRAY * ) createDarray ( 1000, sizeof ( int ) );

	for ( i = 1; i <= num_ctg; i++ )
		{ contig_array[i].flag = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].length + ( unsigned int ) overlaplen >= len_cut )
			{ num_lctg++; }
		else
			{ continue; }

		if ( contig_array[i].flag || contig_array[i].mask || !contig_array[i].downwardConnect || !validConnect ( i, NULL ) )
			{ continue; }

		num5 = num3 = 0;
		ctg = i;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = i;
		contig_array[i].flag = 1;
		bal_ctg = getTwinCtg ( ctg );
		contig_array[bal_ctg].flag = 1;
		len = contig_array[i].length;
		prevCNT = NULL;
		cnt = getNextContig ( ctg, prevCNT, &excep );

		while ( cnt )
		{
			nextCnt = getNextContig ( cnt->contigID, cnt, &excep );

			if ( excep && prevCNT )
				{ fprintf ( stderr, "scaffolding: exception --- prev cnt from %u\n", prevCNT->contigID ); }

			if ( nextCnt && nextCnt->used )
				{ break; }

			setConnectUsed ( ctg, cnt->contigID, 1 );
			* ( int * ) darrayPut ( gap5, num5 - 1 ) = cnt->gapLen;
			ctg = cnt->contigID;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg;
			len += cnt->gapLen + contig_array[ctg].length;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		ctg = getTwinCtg ( i );

		if ( num5 >= 2 )
			{ prevCNT = checkConnect ( getTwinCtg ( * ( unsigned int * ) darrayGet ( scaf5, 1 ) ), ctg ); }
		else
			{ prevCNT = NULL; }

		cnt = getNextContig ( ctg, prevCNT, &excep );

		while ( cnt )
		{
			nextCnt = getNextContig ( cnt->contigID, cnt, &excep );

			if ( excep && prevCNT )
				{ fprintf ( stderr, "scaffolding: exception -- prev cnt from %u\n", prevCNT->contigID ); }

			if ( nextCnt && nextCnt->used )
				{ break; }

			setConnectUsed ( ctg, cnt->contigID, 1 );
			ctg = cnt->contigID;
			len += cnt->gapLen + contig_array[ctg].length;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			* ( int * ) darrayPut ( gap3, num3 ) = cnt->gapLen;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = bal_ctg;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		len += overlaplen;
		sum += len;
		length_array[count++] = len;

		if ( num5 + num3 < 1 )
		{
			fprintf ( stderr, "no scaffold created for contig %d\n", i );
			continue;
		}

		len = prev_ctg = 0;

		for ( j = num3 - 1; j >= 0; j-- )
		{
			if ( prev_ctg )
			{
				num_route = num_trace = 0;
				traceAlongArc ( * ( unsigned int * ) darrayGet ( scaf3, j ), prev_ctg, max_steps
				                , gap - ins_size_var, gap + ins_size_var, 0, 0, &num_route );

				if ( num_route == 1 )
				{
					gap_c++;
				}
			}

			len += contig_array[* ( unsigned int * ) darrayGet ( scaf3, j )].length + * ( int * ) darrayGet ( gap3, j );
			prev_ctg = * ( unsigned int * ) darrayGet ( scaf3, j );
			gap = * ( int * ) darrayGet ( gap3, j ) > 0 ? * ( int * ) darrayGet ( gap3, j ) : 0;
		}

		for ( j = 0; j < num5; j++ )
		{
			if ( prev_ctg )
			{
				num_route = num_trace = 0;
				traceAlongArc ( * ( unsigned int * ) darrayGet ( scaf5, j ), prev_ctg, max_steps
				                , gap - ins_size_var, gap + ins_size_var, 0, 0, &num_route );

				if ( num_route == 1 )
				{
					gap_c++;
				}
			}

			if ( j < num5 - 1 )
			{
				len += contig_array[* ( unsigned int * ) darrayGet ( scaf5, j )].length +
				       * ( int * ) darrayGet ( gap5, j );
				prev_ctg = * ( unsigned int * ) darrayGet ( scaf5, j );
				gap = * ( int * ) darrayGet ( gap5, j ) > 0 ? * ( int * ) darrayGet ( gap5, j ) : 0;
			}
		}
	}

	freeDarray ( scaf3 );
	freeDarray ( scaf5 );
	freeDarray ( gap3 );
	freeDarray ( gap5 );

	if ( count == 0 )
	{
		fprintf ( stderr, "\n\nNo scaffold was constructed.\n\n" );
		free ( ( void * ) length_array );

		for ( j = 0; j < max_n_routes; j++ )
			{ free ( ( void * ) found_routes[j] ); }

		free ( ( void * ) found_routes );
		free ( ( void * ) so_far );
		return;
	}

	fprintf ( stderr, "\nRank %d\n", rank );
	fprintf ( stderr, " Scaffold number                  %d\n", count );
	fprintf ( stderr, " In-scaffold contig number        %u\n", num_lctg / 2 );
	fprintf ( stderr, " Total scaffold length            %lld\n", sum );
	fprintf ( stderr, " Average scaffold length          %lld\n", sum / count );
	fprintf ( stderr, " Filled gap number                %d\n", gap_c );

	//output singleton
	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].length + ( unsigned int ) overlaplen < len_cut || contig_array[i].flag )
			{ continue; }

		length_array[count++] = contig_array[i].length;
		sum += contig_array[i].length;

		if ( isSmallerThanTwin ( i ) )
			{ i++; }
	}

	long int total_len = sum;
	// calculate N50/N90
	qsort ( length_array, count, sizeof ( length_array[0] ), cmp_int );
	N50 = sum * 0.5;
	N90 = sum * 0.9;
	int N50length = 0;
	int N90length = 0;
	sum = flag = 0;

	for ( j = count - 1; j >= 0; j-- )
	{
		sum += length_array[j];

		if ( !flag && sum >= N50 && N50length == 0 )
		{
			N50length = length_array[j];
			flag++;
		}

		if ( sum >= N90 && N90length == 0 )
		{
			N90length = length_array[j];
			break;
		}
	}

	fprintf ( stderr, " Longest scaffold                 %lld\n", length_array[count - 1] );
	fprintf ( stderr, " Scaffold and singleton number    %d\n", count );
	fprintf ( stderr, " Scaffold and singleton length    %lld\n", total_len );
	fprintf ( stderr, " Average length                   %d\n", total_len / count );
	fprintf ( stderr, " N50                              %d\n", N50length );
	fprintf ( stderr, " N90                              %d\n", N90length );
	free ( ( void * ) length_array );

	for ( j = 0; j < max_n_routes; j++ )
		{ free ( ( void * ) found_routes[j] ); }

	free ( ( void * ) found_routes );
	free ( ( void * ) so_far );
}


/*************************************************
Function:
    outputLinks
Description:
    Outputs connection information of current insert size LIB.
Input:
    1. fp:          output file
    2. insertS:     current insert size
Output:
    None.
Return:
    None.
*************************************************/
static void outputLinks ( FILE * fp, int insertS )
{
	unsigned int i, bal_ctg, bal_toCtg;
	CONNECT * cnts, *temp_cnt;

	for ( i = 1; i <= num_ctg; i++ )
	{
		cnts = contig_array[i].downwardConnect;
		bal_ctg = getTwinCtg ( i );

		while ( cnts )
		{
			if ( cnts->weight < 1 )
			{
				cnts = cnts->next;
				continue;
			}

			fprintf ( fp, "%-10d %-10d\t%d\t%d\t%d\n"
			          , i, cnts->contigID, cnts->gapLen, cnts->weight, insertS );
			cnts->weight = 0;
			bal_toCtg = getTwinCtg ( cnts->contigID );
			temp_cnt = getCntBetween ( bal_toCtg, bal_ctg );

			if ( temp_cnt )
				{ temp_cnt->weight = 0; }

			cnts = cnts->next;
		}
	}
}

/*************************************************
 Function:
    PE2Links
 Description:
    Updates connections between contigs based on alignment
    information of paired-end reads.
 Input:
    1. infile:      alignment information file of paired-end reads
 Output:
    None.
 Return:
    None.
 *************************************************/
void PE2Links ( char * infile )
{
	char name[256], *line;
	FILE * fp1;
	FILE * linkF;
	gzFile * fp2;
	int i;
	int flag = 0;
	unsigned int j;
	sprintf ( name, "%s.links", infile );
	boolean filesOK = check_file ( name );

	if ( filesOK )
	{
		fprintf ( stderr, "File %s exists, skip creating the links...\n", name );
		return;
	}

	linkF = ckopen ( name, "w" );

	if ( !pes )
		{ loadPEgrads ( infile ); }

	fprintf ( stderr, "*****************************************************\nStart to load paired-end reads information.\n\n" );

	if ( COMPATIBLE_MODE == 1 )
	{
		sprintf ( name, "%s.readOnContig", infile );
		fp1 = ckopen ( name, "r" );
	}
	else
	{
		sprintf ( name, "%s.readOnContig.gz", infile );
		fp2 = gzopen ( name, "r" );
	}

	lineLen = 1024;
	line = ( char * ) ckalloc ( lineLen * sizeof ( char ) );

	if ( COMPATIBLE_MODE == 1 )
	{
		fgets ( line, lineLen, fp1 );
	}
	else
	{
		gzgets ( fp2, line, lineLen );
	}

	line[0] = '\0';

	for ( i = 0; i < gradsCounter; i++ )
	{
		createCntMemManager();
		createCntLookupTable();
		newCntCounter = 0;

		if ( COMPATIBLE_MODE == 1 )
		{
			flag += connectByPE_grad ( fp1, i, line );
		}
		else
		{
			flag += connectByPE_grad_gz ( fp2, i, line );
		}

		fprintf ( stderr, "%lld new connections.\n\n", newCntCounter / 2 );

		if ( !flag )
		{
			destroyConnectMem();
			deleteCntLookupTable();

			for ( j = 1; j <= num_ctg; j++ )
				{ contig_array[j].downwardConnect = NULL; }

			fprintf ( stderr, "\n" );
			continue;
		}

		flag = 0;
		outputLinks ( linkF, pes[i].insertS );
		destroyConnectMem();
		deleteCntLookupTable();

		for ( j = 1; j <= num_ctg; j++ )
			{ contig_array[j].downwardConnect = NULL; }
	}

	free ( ( void * ) line );

	if ( COMPATIBLE_MODE == 1 )
	{
		fclose ( fp1 );
	}
	else
	{
		gzclose ( fp2 );
	}

	fclose ( linkF );
	fprintf ( stderr, "All paired-end reads information loaded.\n" );
}

/*************************************************
 Function:
    inputLinks
 Description:
    Loads links information of certain insert size.
 Input:
    1. fp:          alignment inofrmation file
    2. insertS:     insert size boundary
    3. line:            the first alignment record of current LIB
 Output:
    1. line:        the first alignment record of next LIB
 Return:
    loaded record number.
 *************************************************/
static int inputLinks ( FILE * fp, int insertS, char * line )
{
	unsigned int ctg, bal_ctg, toCtg, bal_toCtg;
	int gap, wt, ins;
	unsigned int counter = 0, onScafCounter = 0;
	unsigned int maskCounter = 0;
	CONNECT * cnt, *bal_cnt;

	if ( strlen ( line ) )
	{
		sscanf ( line, "%d %d %d %d %d", &ctg, &toCtg, &gap, &wt, &ins );

		if ( ins != insertS )
			{ return counter; }

		if ( 1 )
		{
			bal_ctg = getTwinCtg ( ctg );
			bal_toCtg = getTwinCtg ( toCtg );
			cnt = add1Connect ( ctg, toCtg, gap, wt, 0 );
			bal_cnt = add1Connect ( bal_toCtg, bal_ctg, gap, wt, 0 );

			if ( cnt && insertS > 1000 )
			{
				cnt->newIns = bal_cnt->newIns = 1;
			}

			counter++;

			if ( contig_array[ctg].mask || contig_array[toCtg].mask )
				{ maskCounter++; }

			if ( insertS > 1000 &&
			        contig_array[ctg].from_vt == contig_array[toCtg].from_vt && // on the same scaff
			        contig_array[ctg].indexInScaf < contig_array[toCtg].indexInScaf )
			{
				add1LongPEcov ( ctg, toCtg, wt );
				onScafCounter++;
			}
		}
	}

	while ( fgets ( line, lineLen, fp ) != NULL )
	{
		sscanf ( line, "%d %d %d %d %d", &ctg, &toCtg, &gap, &wt, &ins );

		if ( ins > insertS )
			{ break; }

		if ( insertS > 1000 &&
		        contig_array[ctg].from_vt == contig_array[toCtg].from_vt && // on the same scaff
		        contig_array[ctg].indexInScaf < contig_array[toCtg].indexInScaf )
		{
			add1LongPEcov ( ctg, toCtg, wt );
			onScafCounter++;
		}

		bal_ctg = getTwinCtg ( ctg );
		bal_toCtg = getTwinCtg ( toCtg );
		cnt = add1Connect ( ctg, toCtg, gap, wt, 0 );
		bal_cnt = add1Connect ( bal_toCtg, bal_ctg, gap, wt, 0 );

		if ( cnt && insertS > 1000 )
		{
			cnt->newIns = bal_cnt->newIns = 1;
		}

		counter++;

		if ( contig_array[ctg].mask || contig_array[toCtg].mask )
			{ maskCounter++; }
	}

	fprintf ( stderr, "***************************\nFor insert size: %d\n", insertS );
	fprintf ( stderr, " Total PE links                %d\n", counter );
	fprintf ( stderr, " PE links to masked contigs    %d\n", maskCounter );
	fprintf ( stderr, " On same scaffold PE links     %d\n", onScafCounter );
	return counter;
}

/*************************************************
Function:
    Links2Scaf
Description:
    Constructs scaffolds based on alignment information.
Input:
    1. infile:      prefix of graph
Output:
    None.
Return:
    None.
*************************************************/
void Links2Scaf ( char * infile )
{
	char name[256], *line;
	FILE * fp;
	int i, j = 1, lib_n = 0, cutoff_sum = 0;
	int flag = 0, flag2;
	boolean downS, nonLinear = 0, smallPE = 0, isPrevSmall = 0, markSmall;

	if ( cvg4SNP > 0.001 )
	{
		sprintf ( name, "%s.bubbleInScaff", infile );
		snp_fp = ckopen ( name, "w" );
	}

	cvg4SNP = ( double ) ( cvg4SNP * cvgAvg );

	if ( !pes )
		{ loadPEgrads ( infile ); }

	sprintf ( name, "%s.links", infile );
	fp = ckopen ( name, "r" );
	createCntMemManager();
	createCntLookupTable();
	lineLen = 1024;
	line = ( char * ) ckalloc ( lineLen * sizeof ( char ) );
	fgets ( line, lineLen, fp );
	line[0] = '\0';
	solidArray = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	tempArray = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	scaf3 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	scaf5 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	gap3 = ( DARRAY * ) createDarray ( 1000, sizeof ( int ) );
	gap5 = ( DARRAY * ) createDarray ( 1000, sizeof ( int ) );
	weakPE = 3;
	fprintf ( stderr, "\n" );

	for ( i = 0; i < gradsCounter; i++ )
	{
		if ( MinWeakCut == 0 && i == 0 )
			{ MinWeakCut = pes[i].pair_num_cut; }

		if ( pes[i].insertS < 1000 )
		{
			isPrevSmall = 1;

			if ( MinWeakCut > pes[i].pair_num_cut )
				{ MinWeakCut = pes[i].pair_num_cut; }
		}
		else if ( pes[i].insertS > 1000 && isPrevSmall )
		{
			smallScaf();
			isPrevSmall = 0;
		}

		Insert_size = pes[i].insertS;
		discardCntCounter = 0;
		flag2 = inputLinks ( fp, pes[i].insertS, line );

		if ( flag2 )
		{
			lib_n++;
			cutoff_sum += pes[i].pair_num_cut;
		}

		flag += flag2;

		if ( !flag )
		{
			fprintf ( stderr, "\n" );
			continue;
		}

		if ( i == gradsCounter - 1 || pes[i + 1].rank != pes[i].rank )
		{
			flag = nonLinear = downS = markSmall = 0;

			if ( pes[i].insertS > 1000 && pes[i].rank > 1 )
				{ downS = 1; }

			if ( pes[i].insertS <= 1000 )
				{ smallPE = 1; }

			if ( pes[i].insertS >= 1000 )
			{
				ins_size_var = 50;
				OverlapPercent = 0.05;
			}
			else if ( pes[i].insertS >= 300 )
			{
				ins_size_var = 30;
				OverlapPercent = 0.05;
			}
			else
			{
				ins_size_var = 20;
				OverlapPercent = 0.05;
			}

			if ( pes[i].insertS > 1000 )
				{ weakPE = 5; }

			bySmall = Insert_size > 1000 ? 0 : 1;

			if ( lib_n > 0 )
			{
				weakPE = weakPE < cutoff_sum / lib_n ? cutoff_sum / lib_n : weakPE;
				lib_n = cutoff_sum = 0;
			}

			if ( MinWeakCut > weakPE )
				{ MinWeakCut = weakPE; }

			fprintf ( stderr, "Cutoff of PE links to make a reliable connection: %d\n", weakPE );

			if ( i == gradsCounter - 1 )
				{ nonLinear = 1; }

			if ( Insert_size > 1000 )
			{
				detectBreakScaff();
			}

			ordering ( 1, downS, nonLinear, infile );

			if ( i == gradsCounter - 1 )
			{
				recoverMask();
			}
			else
			{
				scaffold_count ( j, 100 );
				j++;
				fprintf ( stderr, "\n" );
			}

			if ( Insert_size > 1000 && i != gradsCounter - 1 )
			{
				clearNewInsFlag();
			}
		}
	}

	freeDarray ( tempArray );
	freeDarray ( solidArray );
	freeDarray ( scaf3 );
	freeDarray ( scaf5 );
	freeDarray ( gap3 );
	freeDarray ( gap5 );
	free ( ( void * ) line );
	fclose ( fp );

	if ( cvg4SNP > 0.001 )
	{
		fclose ( snp_fp );
	}

	fprintf ( stderr, "\nAll links loaded.\n" );
}

/*************************************************
Function:
    putNodeInArray
Description:
    Puts contig in "ctg4heapArray".
Input:
    1. node:            contig
    2. maxNodes:        maximum allowed contig number in array
    3. dis:         contig's distance to base contig
Output:
    None.
Return:
    1 if contig was already in array or putting operation succeeds.
    0 if array size was larger than cutoff or contig's reverse complement
       format was already in array.
*************************************************/
static boolean putNodeInArray ( unsigned int node, int maxNodes, int dis )
{
	if ( contig_array[node].inSubGraph )
		{ return 1; }

	int index = nodeCounter;

	if ( index > maxNodes )
		{ return 0; }

	if ( contig_array[getTwinCtg ( node )].inSubGraph )
		{ return 0; }

	ctg4heapArray[index].ctgID = node;
	ctg4heapArray[index].dis = dis;
	contig_array[node].inSubGraph = 1;
	ctg4heapArray[index].ds_shut4dheap = 0;
	ctg4heapArray[index].us_shut4dheap = 0;
	ctg4heapArray[index].ds_shut4uheap = 0;
	ctg4heapArray[index].us_shut4uheap = 0;
	return 1;
}

/*************************************************
Function:
    setInGraph
Description:
    Sets contig's status of "inSubGraph".
Input:
    1. flag:        new status.
Output:
    None.
Return:
    None.
*************************************************/
static void setInGraph ( boolean flag )
{
	int i;
	int node;
	nodeCounter = nodeCounter > MaxNodeInSub ? MaxNodeInSub : nodeCounter;

	for ( i = 1; i <= nodeCounter; i++ )
	{
		node = ctg4heapArray[i].ctgID;

		if ( node > 0 )
			{ contig_array[node].inSubGraph = flag; }
	}
}


/*************************************************
Function:
    getIndexInArr
Description:
    Gets contig's index in "ctg4heapArray".
Input:
    1. node:        contig.
Output:
    None.
Return:
    Contig's index if contig was in array.
    0 otherwise.
*************************************************/
static int getIndexInArr ( const unsigned int node )
{
	int i = 1;

	for ( ; i <= nodeCounter; ++i )
	{
		if ( node == ctg4heapArray[i].ctgID )
			{ return i; }
	}

	return 0;
}


/*************************************************
Function:
    dispatch1node
Description:
    Dispatchs a contig to heap according to its distance to base
    contig, then updates related information of heaps.
Input:
    1. dis:         distance to base contig
    2. tempNode:        contig to be dispatched
    3. maxNodes:        maximum allowed contig number in sub-graph
    4. dheap:           heap for downstream contigs of base contig
    5. uheap:           heap for upstream contigs of base contig
    6. DmaxDis:     maximum distance of downstream contig to base contig
    7. UmaxDis:     maximum distance of upstream contig to base
        contig
Output:
    1. dheap:           updated heap for downstream contigs
    2. uheap:           updated heap for upstream contigs
    3. DmaxDis:     updated maximum distance of downstream contig to base contig
    4. UmaxDis:     updated maximum distance of upstream contig to base contig
Return:
    1 if contig was dispatched to "dheap".
    -1 if contig was dispatched to "uheap".
    0 otherwise.
*************************************************/
static boolean dispatch1node ( int dis, unsigned int tempNode, int maxNodes,
                               FibHeap * dheap, FibHeap * uheap, int * DmaxDis, int * UmaxDis )
{
	boolean ret;

	if ( dis >= 0 ) // put it to Dheap
	{
		nodeCounter++;
		ret = putNodeInArray ( tempNode, maxNodes, dis );

		if ( !ret )
			{ return 0; }

		insertNodeIntoHeap ( dheap, dis, nodeCounter );

		if ( dis > *DmaxDis )
			{ *DmaxDis = dis; }

		return 1;
	}
	else           // put it to Uheap
	{
		nodeCounter++;
		ret = putNodeInArray ( tempNode, maxNodes, dis );

		if ( !ret )
			{ return 0; }

		insertNodeIntoHeap ( uheap, -dis, nodeCounter );
		int temp_len = contig_array[tempNode].length;

		if ( -dis + temp_len > *UmaxDis )
			{ *UmaxDis = -dis + contig_array[tempNode].length; }

		return -1;
	}

	return 0;
}

/*************************************************
Function:
    canDheapWait
Description:
    Checks whether current contig is the furthest contig to base
    contig in 'dheap'.
Input:
    1. currNode:        current contig
    2. dis:         current contig's distance to base contig
    3. DmaxDis:     maximum distance of downstream contig to base contig
Output:
    None.
Return:
    0 if current contig was not the furthest contig.
*************************************************/
static boolean canDheapWait ( unsigned int currNode, int dis, int DmaxDis )
{
	if ( dis < DmaxDis )
		{ return 0; }
	else
		{ return 1; }
}

/*************************************************
Function:
    workOnDheap
Description:
    Travers 'dheap' and adds related contigs into 'dheap' or 'uheap'.
Input:
    1. dheap:           heap for downstream contigs of base contig
    2. uheap:           heap for upstream contigs of base contig
    3. Dwait:           indicator of whether all contigs in 'dheap' have been
                            traversed, 1 for yes
    4. Uwait:           indicator of whether all contigs in 'uheap' have been
                         traversed, 1 for yes
    5. DmaxDis:     maximum distance of downstream contig to base
        contig
    6. UmaxDis:     maximum distance of upstream contig to base
        contig
    7. maxNodes:        maximum allowed contig number in sub-graph
Output:
    1. dheap:           updated heap for downstream contigs
    2. uheap:           updated heap for upstream contigs
    3. Dwait:           indicator of whether all contigs in 'dheap' have been
                            traversed, 1 for yes
    4. Uwait:           indicator of whether all contigs in 'uheap' have been
                            traversed, 1 for yes
    5. DmaxDis:     updated maximum distance of downstream contig to
                        base contig
    6. UmaxDis:     updated maximum distance of upstream contig to
                        base contig
Return:
    0 if operation of putting contig into sub-graph failed.
*************************************************/
static boolean workOnDheap ( FibHeap * dheap, FibHeap * uheap, boolean * Dwait, boolean * Uwait,
                             int * DmaxDis, int * UmaxDis, int maxNodes )
{
	if ( *Dwait )
		{ return 1; }

	unsigned int currNode, twin, tempNode;
	CTGinHEAP * ctgInHeap;
	int indexInArray;
	CONNECT * us_cnt, *ds_cnt;
	int dis0, dis;
	boolean ret, isEmpty;

	while ( ( indexInArray = removeNextNodeFromHeap ( dheap ) ) != 0 )
	{
		ctgInHeap = &ctg4heapArray[indexInArray];
		currNode = ctgInHeap->ctgID;
		dis0 = ctgInHeap->dis;
		isEmpty = IsHeapEmpty ( dheap );
		twin = getTwinCtg ( currNode );
		us_cnt = ctgInHeap->us_shut4dheap ? NULL : contig_array[twin].downwardConnect;

		while ( us_cnt )
		{
			if ( us_cnt->deleted || us_cnt->mask ||
			        contig_array[getTwinCtg ( us_cnt->contigID )].inSubGraph )
			{
				us_cnt = us_cnt->next;
				continue;
			}

			tempNode = getTwinCtg ( us_cnt->contigID );

			if ( contig_array[tempNode].inSubGraph )
			{
				us_cnt = us_cnt->next;
				continue;
			}

			dis = dis0 - us_cnt->gapLen - ( int ) contig_array[twin].length;
			ret = dispatch1node ( dis, tempNode, maxNodes, dheap, uheap, DmaxDis, UmaxDis );

			if ( ret == 0 )
				{ return 0; }
			else if ( ret < 0 )
				{ *Uwait = 0; }

			us_cnt = us_cnt->next;
		}

		if ( nodeCounter > 1 && isEmpty )
		{
			*Dwait = canDheapWait ( currNode, dis0, *DmaxDis );

			if ( *Dwait )
			{
				isEmpty = IsHeapEmpty ( dheap );
				insertNodeIntoHeap ( dheap, dis0, indexInArray );
				ctg4heapArray[indexInArray].us_shut4dheap = 1;

				if ( isEmpty )
					{ return 1; }
				else
					{ continue; }
			}
		}

		ds_cnt = ctgInHeap->ds_shut4dheap ? NULL : contig_array[currNode].downwardConnect;

		while ( ds_cnt )
		{
			if ( ds_cnt->deleted || ds_cnt->mask || contig_array[ds_cnt->contigID].inSubGraph )
			{
				ds_cnt = ds_cnt->next;
				continue;
			}

			tempNode = ds_cnt->contigID;
			dis = dis0 + ds_cnt->gapLen + ( int ) contig_array[tempNode].length;
			ret = dispatch1node ( dis, tempNode, maxNodes, dheap, uheap, DmaxDis, UmaxDis );

			if ( ret == 0 )
				{ return 0; }
			else if ( ret < 0 )
				{ *Uwait = 0; }

			ds_cnt = ds_cnt->next;
		}  // for each downstream connections
	}  // for each node comes off the heap

	*Dwait = 1;
	return 1;
}

/*************************************************
Function:
    canUheapWait
Description:
    Checks whether current contig is the furthest contig to base
    contig in "uheap".
Input:
    1. currNode:        current contig
    2. dis:         current contig's distance to base contig
    3. DmaxDis:     maximum distance of downstream contig to base contig
Output:
    None.
Return:
    0 if current contig was not the furthest contig.
*************************************************/
static boolean canUheapWait ( unsigned int currNode, int dis, int UmaxDis )
{
	int temp_len = contig_array[currNode].length;

	if ( -dis + temp_len < UmaxDis )
		{ return 0; }
	else
		{ return 1; }
}

/*************************************************
Function:
    workOnUheap
Description:
    Travers 'uheap' and adds related contigs into "dheap" or "uheap".
Input:
    1. dheap:           heap for downstream contigs of base contig
    2. uheap:           heap for upstream contigs of base contig
    3. Dwait:           indicator of whether all contigs in 'dheap' have been
                            traversed, 1 for yes
    4. Uwait:           indicator of whether all contigs in 'uheap' have been
                            traversed, 1 for yes
    5. DmaxDis:     maximum distance of downstream contig to base contig
    6. UmaxDis:     maximum distance of upstream contig to base contig
    7. maxNodes:        maximum allowed contig number in sub-graph
Output:
    1. dheap:           updated heap for downstream contigs
    2. uheap:           updated heap for upstream contigs
    3. Dwait:           indicator of whether all contigs in "dheap" have been
                        traversed, 1 for yes
    4. Uwait:           indicator of whether all contigs in "uheap" have been
                        traversed, 1 for yes
    5. DmaxDis:     updated maximum distance of downstream contig to
                        base contig
    6. UmaxDis:     updated maximum distance of upstream contig to
                        base contig
Return:
    0 if operation of putting contig into sub-graph failed.
*************************************************/
static boolean workOnUheap ( FibHeap * dheap, FibHeap * uheap, boolean * Dwait, boolean * Uwait,
                             int * DmaxDis, int * UmaxDis, int maxNodes )
{
	if ( *Uwait )
		{ return 1; }

	unsigned int currNode, twin, tempNode;
	CTGinHEAP * ctgInHeap;
	int indexInArray;
	CONNECT * us_cnt, *ds_cnt;
	int dis0, dis;
	boolean ret, isEmpty;

	while ( ( indexInArray = removeNextNodeFromHeap ( uheap ) ) != 0 )
	{
		ctgInHeap = &ctg4heapArray[indexInArray];
		currNode = ctgInHeap->ctgID;
		dis0 = ctgInHeap->dis;
		isEmpty = IsHeapEmpty ( uheap );
		ds_cnt = ctgInHeap->ds_shut4uheap ? NULL : contig_array[currNode].downwardConnect;

		while ( ds_cnt )
		{
			if ( ds_cnt->deleted || ds_cnt->mask || contig_array[ds_cnt->contigID].inSubGraph )
			{
				ds_cnt = ds_cnt->next;
				continue;
			}

			tempNode = ds_cnt->contigID;
			dis = dis0 + ds_cnt->gapLen + contig_array[tempNode].length;
			ret = dispatch1node ( dis, tempNode, maxNodes, dheap, uheap, DmaxDis, UmaxDis );

			if ( ret == 0 )
				{ return 0; }
			else if ( ret > 0 )
				{ *Dwait = 0; }

			ds_cnt = ds_cnt->next;
		}  // for each downstream connections

		if ( nodeCounter > 1 && isEmpty )
		{
			*Uwait = canUheapWait ( currNode, dis0, *UmaxDis );

			if ( *Uwait )
			{
				isEmpty = IsHeapEmpty ( uheap );
				insertNodeIntoHeap ( uheap, dis0, indexInArray );
				ctg4heapArray[indexInArray].ds_shut4uheap = 1;

				if ( isEmpty )
					{ return 1; }
				else
					{ continue; }
			}
		}

		twin = getTwinCtg ( currNode );
		us_cnt = ctgInHeap->us_shut4uheap ? NULL : contig_array[twin].downwardConnect;

		while ( us_cnt )
		{
			if ( us_cnt->deleted || us_cnt->mask ||
			        contig_array[getTwinCtg ( us_cnt->contigID )].inSubGraph )
			{
				us_cnt = us_cnt->next;
				continue;
			}

			tempNode = getTwinCtg ( us_cnt->contigID );

			if ( contig_array[tempNode].inSubGraph )
			{
				us_cnt = us_cnt->next;
				continue;
			}

			dis = dis0 - us_cnt->gapLen - contig_array[twin].length;
			ret = dispatch1node ( dis, tempNode, maxNodes, dheap, uheap, DmaxDis, UmaxDis );

			if ( ret == 0 )
				{ return 0; }
			else if ( ret > 0 )
				{ *Dwait = 0; }

			us_cnt = us_cnt->next;
		}
	}  // for each node comes off the heap

	*Uwait = 1;
	return 1;
}

/*************************************************
Function:
    pickUpGeneralSubgraph
Description:
    Starting from a base contig, gathers related contigs to form a
    sub-graph.
Input:
    1. node1:           base contig
    2. maxNodes:        maximum allowed contig number in sub-graph
Output:
    None.
Return:
    0 if sub-graph was not found.
*************************************************/
static boolean pickUpGeneralSubgraph ( unsigned int node1, int maxNodes )
{
	FibHeap * Uheap = newFibHeap(); // heap for upstream contigs to node1
	FibHeap * Dheap = newFibHeap();
	int UmaxDis;   // max distance upstream to node1
	int DmaxDis;
	boolean Uwait;  // wait signal for Uheap
	boolean Dwait;
	int dis;
	boolean ret;
	//initiate: node1 is put to array once, and to both Dheap and Uheap
	dis = 0;
	nodeCounter = 1;
	putNodeInArray ( node1, maxNodes, dis );
	insertNodeIntoHeap ( Dheap, dis, nodeCounter );
	ctg4heapArray[nodeCounter].us_shut4dheap = 1;
	Dwait = 0;
	DmaxDis = 0;
	insertNodeIntoHeap ( Uheap, dis, nodeCounter );
	ctg4heapArray[nodeCounter].ds_shut4uheap = 1;
	Uwait = 1;
	UmaxDis = contig_array[node1].length;

	while ( 1 )
	{
		ret = workOnDheap ( Dheap, Uheap, &Dwait, &Uwait, &DmaxDis, &UmaxDis, maxNodes );

		if ( !ret )
		{
			setInGraph ( 0 );
			destroyHeap ( Dheap );
			destroyHeap ( Uheap );
			return 0;
		}

		ret = workOnUheap ( Dheap, Uheap, &Dwait, &Uwait, &DmaxDis, &UmaxDis, maxNodes );

		if ( !ret )
		{
			setInGraph ( 0 );
			destroyHeap ( Dheap );
			destroyHeap ( Uheap );
			return 0;
		}

		if ( Uwait && Dwait )
		{
			destroyHeap ( Dheap );
			destroyHeap ( Uheap );
			return 1;
		}
	}
}

/*************************************************
Function:
    cmp_ctg
Description:
    Compares two contigs according to their distances to base contig.
Input:
    1. a:       pointer to the first contig in heap
    2. b:       pointer to the second contig in heap
Output:
    None.
Return:
    1 if the first contig was further than the second contig.
    -1 if the second contig was further than the first contig.
    0 if the distances were equal.
*************************************************/
static int cmp_ctg ( const void * a, const void * b )
{
	CTGinHEAP * A, *B;
	A = ( CTGinHEAP * ) a;
	B = ( CTGinHEAP * ) b;

	if ( A->dis > B->dis )
		{ return 1; }
	else if ( A->dis == B->dis )
		{ return 0; }
	else
		{ return -1; }
}

/*************************************************
Function:
    checkEligible
Description:
    Checks if the sub-graph accords with the follwing criterions:
    1. The reversed complement of the first contig has neither
        downstream connection to contigs in sub-graph, nor more
        than one upstream connections in scaffold.
    2. The last contig has neither downstream connection to reversed
        complement of any contig in sub-graph, nor more than one
        upstream connections in scaffold.
    3. Except the first and the last contigs, none of other contigs in
        sub-graph has connection to contig which is out of sub-graph.
Input:
    None.
Output:
    None.
Return:
    1 if the sub-graph accorded with all above criterions.
*************************************************/
static boolean checkEligible()
{
	unsigned int firstNode = ctg4heapArray[1].ctgID;
	unsigned int twin;
	int i;
	boolean flag = 0;
	//check if the first node has incoming link from twin of any node in subgraph
	// or it has multi outgoing links bound to incoming links
	twin = getTwinCtg ( firstNode );
	CONNECT * ite_cnt = contig_array[twin].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->deleted || ite_cnt->mask )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		if ( contig_array[ite_cnt->contigID].inSubGraph )
		{
			return 0;
		}

		if ( ite_cnt->prevInScaf )
		{
			if ( flag )
				{ return 0; }

			flag = 1;
		}

		ite_cnt = ite_cnt->next;
	}

	//check if the last node has outgoing link to twin of any node in subgraph
	// or it has multi outgoing links bound to incoming links
	unsigned int lastNode = ctg4heapArray[nodeCounter].ctgID;
	ite_cnt = contig_array[lastNode].downwardConnect;
	flag = 0;

	while ( ite_cnt )
	{
		if ( ite_cnt->deleted || ite_cnt->mask )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		twin = getTwinCtg ( ite_cnt->contigID );

		if ( contig_array[twin].inSubGraph )
		{
			return 0;
		}

		if ( ite_cnt->prevInScaf )
		{
			if ( flag )
				{ return 0; }

			flag = 1;
		}

		ite_cnt = ite_cnt->next;
	}

	//check if any node has outgoing link to node outside the subgraph
	for ( i = 1; i < nodeCounter; i++ )
	{
		ite_cnt = contig_array[ctg4heapArray[i].ctgID].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->deleted || ite_cnt->mask )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			if ( !contig_array[ite_cnt->contigID].inSubGraph )
			{
				return 0;
			}

			ite_cnt = ite_cnt->next;
		}
	}

	//check if any node has incoming link from node outside the subgraph
	for ( i = 2; i <= nodeCounter; i++ )
	{
		twin = getTwinCtg ( ctg4heapArray[i].ctgID );
		ite_cnt = contig_array[twin].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->deleted || ite_cnt->mask )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			if ( !contig_array[getTwinCtg ( ite_cnt->contigID )].inSubGraph )
			{
				return 0;
			}

			ite_cnt = ite_cnt->next;
		}
	}

	return 1;
}

/*************************************************
Function:
    arrayvalue
Description:
    Copies one contig's values to another in array.
Input:
    1. init_array:          pointer to target contig in array
    2. value_array:     pointer to source contig in array
Output:
    None.
Return:
    None.
*************************************************/
static void arrayvalue ( CTGinHEAP * init_array, CTGinHEAP * value_array )
{
	init_array->ctgID = value_array->ctgID;
	init_array->dis = value_array->dis;
	init_array->ds_shut4dheap = value_array->ds_shut4dheap;
	init_array->ds_shut4uheap = value_array->ds_shut4uheap;
	init_array->us_shut4dheap = value_array->us_shut4dheap;
	init_array->us_shut4uheap = value_array->us_shut4uheap;
}

/*************************************************
Function:
    arrayexchange
Description:
    Exchanges two contigs' values.
Input:
    1. from_id:     the first contig in array
    2. to_id:           the second contig in array
Output:
    None.
Return:
    None.
*************************************************/
static void arrayexchange ( unsigned int from_id, unsigned int to_id )
{
	CTGinHEAP tmp_h;
	arrayvalue ( &tmp_h, & ( ctg4heapArray[from_id] ) );
	arrayvalue ( & ( ctg4heapArray[from_id] ), & ( ctg4heapArray[to_id] ) );
	arrayvalue ( & ( ctg4heapArray[to_id] ), &tmp_h );
}

static void deletearray ( unsigned int id )
{
	int i;

	for ( i = 1; i < nodeCounter; i++ )
	{
		if ( i >= id )
		{
			arrayvalue ( & ( ctg4heapArray[i] ), & ( ctg4heapArray[i + 1] ) );
		}
	}

	nodeCounter--;
}

int getnextInScafCtg ( int id, int mask, int flag )
{
	CONNECT * tmp_cn = contig_array[id].downwardConnect;
	int currId = 0;

	while ( tmp_cn )
	{
		if ( tmp_cn->prevInScaf )
		{
			currId = tmp_cn->contigID;

			if ( mask != 0 && tmp_cn->contigID != mask )
				{ break; }
		}

		tmp_cn = tmp_cn->next;
	}

	if ( mask == currId && mask != 0 )
		{ currId = 0; }

	if ( flag == 0 && currId != 0 )
		{ currId = getTwinCtg ( currId ); }

	return currId;
}

void delete_PrevNext ( int i, int flag )
{
	int id = ctg4heapArray[i].ctgID;
	int pid = getTwinCtg ( id );
	CONNECT * ite_cnt = contig_array[id].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->mask || ite_cnt->deleted || !contig_array[ite_cnt->contigID].inSubGraph )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		if ( flag == 1 )
			{ setNextInScaf ( ite_cnt, NULL ); }

		if ( flag == 0 )
			{ setPrevInScaf ( ite_cnt, 0 ); }

		ite_cnt = ite_cnt->next;
	}

	ite_cnt = contig_array[pid].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->mask || ite_cnt->deleted || !contig_array[ite_cnt->contigID].inSubGraph )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		if ( flag == 0 )
			{ setNextInScaf ( ite_cnt, NULL ); }

		if ( flag == 1 )
			{ setPrevInScaf ( ite_cnt, 0 ); }

		ite_cnt = ite_cnt->next;
	}
}

/*************************************************
Function:
    getCntNodes
Description:
    Gets contig's downstream connection provided by short insert
    size paired-end reads.
Input:
    1. node:            contig
Output:
    1. nodeArray:       array to store downstream contigs
    2. gapArray:        array to store distance from current contig to
                        downstream contigs
Return:
    Found downstream contig number.
*************************************************/
static int getCntNodes ( unsigned int node, unsigned int * nodeArray, unsigned int * gapArray )
{
	int count = 0;
	CONNECT * cnt = contig_array[node].downwardConnect;

	while ( cnt )
	{
		if ( 0 == bySmall && cnt->weight < 3 && !cnt->smallIns && !cnt->bySmall )
		{
			cnt = cnt->next;
			continue;
		}

		nodeArray[count] = cnt->contigID;
		gapArray[count++] = cnt->gapLen;

		if ( count == MaxCntNode )
		{
			break;
		}

		cnt = cnt->next;
	}

	return count;
}

/*************************************************
Function:
    calGapLen
Description:
    Calculates two contigs' distance through the common contig
    to which both contigs connected.
Input:
    1. cntCounter:          contig numbet in the "cntNodeArr"
    2. cntNodeArr:      the array of third-party contigs
    3. cntGapArr:           the array of distances to third-party contigs
    4. node1:               the first contig
    5. node2:               the second contig
    6. tmpNode:         the contig whose upstream/downstream contigs
                            array was not selected as "cntNodeArr"
Output:
    1. cntCounter:          common contig number of these two contigs
Return:
    Accumulated distance between these two contigs.
*************************************************/
static int calGapLen ( int * cntCounter, unsigned int * cntNodeArr, int * cntGapArr,
                       unsigned int node1, unsigned int node2, unsigned int tmpNode )
{
	int i = 0, gapLen = 0, count = 0;
	unsigned int target_node;
	CONNECT * cnt;
	int len = contig_array[node2].length;

	for ( ; i < *cntCounter; ++i )
	{
		cnt = getCntBetween ( tmpNode, cntNodeArr[i] );

		if ( cnt && ( cnt->weight >= 3 || bySmall || cnt->smallIns || cnt->bySmall ) )
		{
			if ( tmpNode == node1 )
			{
				gapLen += cnt->gapLen - cntGapArr[i] - len;
			}
			else
			{
				gapLen += cntGapArr[i] - cnt->gapLen - len;
			}

			++count;
		}
	}

	*cntCounter = count;
	return gapLen;
}

/*************************************************
Function:
    arrangeNodes_general
Description:
    Arranges contigs in sub-graph and updates related relation.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void arrangeNodes_general()
{
	int i, j, gap, adjustedGap;
	CONNECT * ite_cnt, *temp_cnt, *bal_cnt, *prev_cnt, *next_cnt, *dh_cnt, *three_cnt;
	unsigned int node1, node2, tmpNode;
	unsigned int bal_nd1, bal_nd2;
	unsigned int pre_node, bal_pre_node, next_node, bal_next_node;  // pre/next node that connected to first/last node if there is
	unsigned int first_node, bal_first_node, last_node, bal_last_node;
	unsigned int affected_node1, bal_affected_node1, affected_node2, bal_affected_node2;
	unsigned int exchangeNode1 = 0, exchangeNode2 = 0;
	int cntCounter, comCount;
	int tmp_dis;

	//delete original connections
	for ( i = 1; i <= nodeCounter; i++ )
	{
		node1 = ctg4heapArray[i].ctgID;
		ite_cnt = contig_array[node1].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->mask || ite_cnt->deleted || !contig_array[ite_cnt->contigID].inSubGraph )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			ite_cnt->deleted = 1;
			setNextInScaf ( ite_cnt, NULL );
			setPrevInScaf ( ite_cnt, 0 );
			ite_cnt = ite_cnt->next;
		}

		bal_nd1 = getTwinCtg ( node1 );
		ite_cnt = contig_array[bal_nd1].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->mask || ite_cnt->deleted || !contig_array[getTwinCtg ( ite_cnt->contigID )].inSubGraph )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			ite_cnt->deleted = 1;
			setNextInScaf ( ite_cnt, NULL );
			setPrevInScaf ( ite_cnt, 0 );
			ite_cnt = ite_cnt->next;
		}
	}

	CONNECT * first_cnt = NULL, *last_cnt = NULL, *tmp_cnt;
	CONNECT * affected_cnt = NULL, *bal_affected_cnt = NULL;    //connections connected to pre_cnt and next_cnt
	pre_node = bal_pre_node = next_node = bal_next_node = 0;
	first_node = ctg4heapArray[1].ctgID;
	bal_first_node = getTwinCtg ( first_node );
	ite_cnt = contig_array[bal_first_node].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->deleted || ite_cnt->mask )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		if ( ite_cnt->prevInScaf )
		{
			if ( !first_cnt )
			{
				first_cnt = ite_cnt;
			}

			bal_pre_node = ite_cnt->contigID;
			pre_node = getTwinCtg ( bal_pre_node );
			tmp_cnt = getCntBetween ( pre_node, first_node );
		}

		ite_cnt = ite_cnt->next;
	}

	last_node = ctg4heapArray[nodeCounter].ctgID;
	bal_last_node = getTwinCtg ( last_node );
	ite_cnt = contig_array[last_node].downwardConnect;

	while ( ite_cnt )
	{
		if ( ite_cnt->deleted || ite_cnt->mask )
		{
			ite_cnt = ite_cnt->next;
			continue;
		}

		if ( ite_cnt->prevInScaf )
		{
			if ( !last_cnt )
			{
				last_cnt = ite_cnt;
			}

			next_node = ite_cnt->contigID;
			bal_next_node = getTwinCtg ( next_node );
			tmp_cnt = getCntBetween ( bal_next_node, bal_last_node );
		}

		ite_cnt = ite_cnt->next;
	}

	prev_cnt = next_cnt = NULL;

	for ( i = 1; i < nodeCounter; ++i )
	{
		node1 = ctg4heapArray[i].ctgID;
		node2 = ctg4heapArray[i + 1].ctgID;
		bal_nd1 = getTwinCtg ( node1 );
		bal_nd2 = getTwinCtg ( node2 );
		gap = ctg4heapArray[i + 1].dis - ctg4heapArray[i].dis
		      - contig_array[node2].length;
		temp_cnt = getCntBetween ( node1, node2 );
		dh_cnt = getCntBetween ( node2, node1 );

		if ( i >= 2 )
			{ three_cnt = getCntBetween ( ctg4heapArray[i - 1].ctgID, node2 ); }

		if ( dh_cnt )
		{
			tmp_dis = ( int ) contig_array[node1].length + ( int ) contig_array[node2].length + gap + dh_cnt->gapLen;
		}
		else
		{
			tmp_dis = -1;
		}

		if ( temp_cnt && ( bySmall || temp_cnt->bySmall || temp_cnt->smallIns || !dh_cnt || !dh_cnt->bySmall || !dh_cnt->smallIns ) )
		{
			temp_cnt->deleted = 0;
			temp_cnt->mask = 0;
			bal_cnt = getCntBetween ( bal_nd2, bal_nd1 );
			bal_cnt->deleted = 0;
			bal_cnt->mask = 0;
		}
		else if ( dh_cnt && ( ( dh_cnt->bySmall || dh_cnt->smallIns || bySmall )
		                      || ( ( -gap > ( int ) contig_array[node1].length || -gap > ( int ) contig_array[node2].length )
		                           && tmp_dis > 0 && tmp_dis < 500 && dh_cnt->weight > 3 ) ) )
		{
			dh_cnt->deleted = 0;
			dh_cnt->mask = 0;
			dh_cnt = getCntBetween ( bal_nd1, bal_nd2 );
			dh_cnt->deleted = 0;
			dh_cnt->mask = 0;
			arrayexchange ( i, i + 1 );

			if ( i == 1 )
			{
				i = 0;
				continue;
			}

			if ( i == 2 )
			{
				prev_cnt->deleted = 1;
				prev_cnt = NULL;
				next_cnt->deleted = 1;
				next_cnt = NULL;
				i = 0;
				continue;
			}

			bal_affected_node2 = next_cnt->contigID;
			affected_node2 = getTwinCtg ( bal_affected_node2 );
			bal_affected_cnt = next_cnt->nextInScaf;
			bal_affected_node1 = bal_affected_cnt->contigID;
			affected_node1 = getTwinCtg ( bal_affected_node1 );
			affected_cnt = getCntBetween ( affected_node1, affected_node2 );

			if ( !affected_cnt )
			{
				fprintf ( stderr, "affected cnt between %u(%u) and %u(%u) doesn't exists!\n", affected_node1, bal_affected_node1, affected_node2, bal_affected_node2 );
				exit ( 1 );
			}

			setNextInScaf ( affected_cnt, NULL );
			setPrevInScaf ( bal_affected_cnt, 0 );
			setNextInScaf ( next_cnt, NULL );
			setPrevInScaf ( prev_cnt, 0 );
			prev_cnt->deleted = 1;
			prev_cnt = NULL;
			next_cnt->deleted = 1;
			next_cnt = NULL;
			i -= 3;
			continue;
		}
		else
		{
			if ( ( bySmall > 0 && gap < 0 )
			        || ( -gap > ( int ) contig_array[node1].length || -gap > ( int ) contig_array[node2].length )
			        && ( i != nodeCounter - 1 ) )
			{
				adjustedGap = comCount = 0;
				uCntCounter1 = getCntNodes ( getTwinCtg ( node1 ), uCntNodeArr1, uCntGapArr1 );
				dCntCounter1 = getCntNodes ( node1, dCntNodeArr1, dCntGapArr1 );
				uCntCounter2 = getCntNodes ( getTwinCtg ( node2 ), uCntNodeArr2, uCntGapArr2 );
				dCntCounter2 = getCntNodes ( node2, dCntNodeArr2, dCntGapArr2 );

				if ( uCntCounter1 < uCntCounter2 )
				{
					tmpNode = getTwinCtg ( node2 );
					cntCounter = uCntCounter1;
					cntNodeArr = &uCntNodeArr1[0];
					cntGapArr  = &uCntGapArr1[0];
				}
				else
				{
					tmpNode = getTwinCtg ( node1 );
					cntCounter = uCntCounter2;
					cntNodeArr = &uCntNodeArr2[0];
					cntGapArr  = &uCntGapArr2[0];
				}

				adjustedGap += calGapLen ( &cntCounter, cntNodeArr, cntGapArr, getTwinCtg ( node2 ), getTwinCtg ( node1 ), tmpNode );
				comCount += cntCounter;

				if ( dCntCounter1 < dCntCounter2 )
				{
					tmpNode = node2;
					cntCounter = dCntCounter1;
					cntNodeArr = &dCntNodeArr1[0];
					cntGapArr  = &dCntGapArr1[0];
				}
				else
				{
					tmpNode = node1;
					cntCounter = dCntCounter2;
					cntNodeArr = &dCntNodeArr2[0];
					cntGapArr  = &dCntGapArr2[0];
				}

				adjustedGap += calGapLen ( &cntCounter, cntNodeArr, cntGapArr, node1, node2, tmpNode );
				comCount += cntCounter;

				if ( comCount > 0 )
				{
					gap = adjustedGap / comCount;
				}
			}

			if ( ( -gap > ( int ) contig_array[node1].length || -gap > ( int ) contig_array[node2].length )
			        && ( i != nodeCounter - 1 ) && ( ( exchangeNode1 == 0 && exchangeNode2 == 0 )
			                || ( exchangeNode1 != ctg4heapArray[i + 1].ctgID && exchangeNode2 != ctg4heapArray[i].ctgID ) ) )
			{
				exchangeNode1 = ctg4heapArray[i].ctgID;
				exchangeNode2 = ctg4heapArray[i + 1].ctgID;
				arrayexchange ( i, i + 1 );

				if ( i == 1 )
				{
					i--;
					continue;
				}

				if ( i == 2 )
				{
					prev_cnt->deleted = 1;
					prev_cnt = NULL;
					next_cnt->deleted = 1;
					next_cnt = NULL;
					i = 0;
					continue;
				}

				bal_affected_node2 = next_cnt->contigID;
				affected_node2 = getTwinCtg ( bal_affected_node2 );
				bal_affected_cnt = next_cnt->nextInScaf;
				bal_affected_node1 = bal_affected_cnt->contigID;
				affected_node1 = getTwinCtg ( bal_affected_node1 );
				affected_cnt = getCntBetween ( affected_node1, affected_node2 );

				if ( !affected_cnt )
				{
					fprintf ( stderr, "affected cnt between %u(%u) and %u(%u) doesn't exists!\n", affected_node1, bal_affected_node1, affected_node2, bal_affected_node2 );
					exit ( 1 );
				}

				setNextInScaf ( affected_cnt, NULL );
				setPrevInScaf ( bal_affected_cnt, 0 );
				setNextInScaf ( next_cnt, NULL );
				setPrevInScaf ( prev_cnt, 0 );
				prev_cnt->deleted = 1;
				prev_cnt = NULL;
				next_cnt->deleted = 1;
				next_cnt = NULL;
				i -= 3;
				continue;
			}

			temp_cnt = allocateCN ( node2, gap );

			if ( cntLookupTable )
				{ putCnt2LookupTable ( node1, temp_cnt ); }

			temp_cnt->weight = 0;
			temp_cnt->next = contig_array[node1].downwardConnect;
			contig_array[node1].downwardConnect = temp_cnt;
			bal_cnt = allocateCN ( bal_nd1, gap );

			if ( cntLookupTable )
				{ putCnt2LookupTable ( bal_nd2, bal_cnt ); }

			bal_cnt->weight = 0;
			bal_cnt->next = contig_array[bal_nd2].downwardConnect;
			contig_array[bal_nd2].downwardConnect = bal_cnt;
		}

		if ( prev_cnt )
		{
			setNextInScaf ( prev_cnt, temp_cnt );
			setPrevInScaf ( temp_cnt, 1 );
		}

		if ( next_cnt )
		{
			setNextInScaf ( bal_cnt, next_cnt );
			setPrevInScaf ( next_cnt, 1 );
		}

		prev_cnt = temp_cnt;
		next_cnt = bal_cnt;
	}

	if ( first_cnt )
	{
		if ( ctg4heapArray[1].ctgID == first_node )
		{
			bal_nd1 = first_cnt->contigID;
			node1 = getTwinCtg ( bal_nd1 );
			node2 = first_node;
			temp_cnt = checkConnect ( node1, node2 );
			bal_cnt = first_cnt;
			next_cnt = checkConnect ( ctg4heapArray[1].ctgID, ctg4heapArray[2].ctgID );
			prev_cnt = checkConnect ( getTwinCtg ( ctg4heapArray[2].ctgID ), getTwinCtg ( ctg4heapArray[1].ctgID ) );

			if ( temp_cnt )
			{
				setNextInScaf ( temp_cnt, next_cnt );
				setPrevInScaf ( temp_cnt->nextInScaf, 0 );
				setPrevInScaf ( next_cnt, 1 );
				setNextInScaf ( prev_cnt, bal_cnt );
			}
		}
		else
		{
			bal_pre_node = first_cnt->contigID;
			pre_node = getTwinCtg ( bal_pre_node );
			j = 1;
			node1 = ctg4heapArray[j].ctgID;
			node2 = ctg4heapArray[j + 1].ctgID;
			ite_cnt = getCntBetween ( pre_node, node1 );
			bal_cnt = getCntBetween ( getTwinCtg ( node1 ), bal_pre_node );

			while ( !ite_cnt && node2 != first_node )
			{
				tmp_cnt = getCntBetween ( node1, node2 );
				bal_cnt = getCntBetween ( getTwinCtg ( node2 ), getTwinCtg ( node1 ) );
				setNextInScaf ( tmp_cnt, NULL );
				setPrevInScaf ( tmp_cnt, 0 );
				tmp_cnt->deleted = 1;
				setNextInScaf ( bal_cnt, NULL );
				setPrevInScaf ( bal_cnt, 0 );
				bal_cnt->deleted = 1;
				++j;
				node1 = ctg4heapArray[j].ctgID;
				node2 = ctg4heapArray[j + 1].ctgID;
				ite_cnt = getCntBetween ( pre_node, node1 );
				bal_cnt = getCntBetween ( getTwinCtg ( node1 ), bal_pre_node );
			}

			if ( !ite_cnt )
			{
				tmp_cnt = getCntBetween ( node1, first_node );
				gap = first_cnt->gapLen - tmp_cnt->gapLen - contig_array[node1].length;
				ite_cnt = allocateCN ( node1, gap );
				ite_cnt->weight = 0;

				if ( cntLookupTable )
				{
					putCnt2LookupTable ( pre_node, ite_cnt );
				}

				ite_cnt->next = contig_array[pre_node].downwardConnect;
				contig_array[pre_node].downwardConnect = ite_cnt;
				bal_cnt = allocateCN ( bal_pre_node, gap );
				bal_cnt->weight = 0;

				if ( cntLookupTable )
				{
					putCnt2LookupTable ( getTwinCtg ( node1 ), bal_cnt );
				}

				bal_cnt->next = contig_array[getTwinCtg ( node1 )].downwardConnect;
				contig_array[getTwinCtg ( node1 )].downwardConnect = bal_cnt;
			}

			ite_cnt->deleted = 0;
			ite_cnt->mask = 0;
			bal_cnt->deleted = 0;
			bal_cnt->mask = 0;
			tmp_cnt = getCntBetween ( node1, node2 );
			setNextInScaf ( ite_cnt, tmp_cnt );
			setPrevInScaf ( tmp_cnt, 1 );
			tmp_cnt = getCntBetween ( getTwinCtg ( node2 ), getTwinCtg ( node1 ) );
			setNextInScaf ( tmp_cnt, bal_cnt );
			setPrevInScaf ( bal_cnt, 1 );

			if ( first_cnt->nextInScaf )
			{
				setNextInScaf ( bal_cnt, first_cnt->nextInScaf );
				bal_affected_node1 = first_cnt->nextInScaf->contigID;
				affected_node1 = getTwinCtg ( bal_affected_node1 );
				affected_cnt = getCntBetween ( affected_node1, pre_node );
				setNextInScaf ( affected_cnt, ite_cnt );
				setPrevInScaf ( ite_cnt, 1 );
			}

			setNextInScaf ( first_cnt, NULL );
			setPrevInScaf ( first_cnt, 0 );
			first_cnt->deleted = 1;
			first_cnt->mask = 1;
			tmp_cnt = getCntBetween ( pre_node, first_node );
			setNextInScaf ( tmp_cnt, NULL );
			setPrevInScaf ( tmp_cnt, 0 );
			tmp_cnt->deleted = 1;
			tmp_cnt->mask = 1;
		}
	}

	if ( last_cnt )
	{
		node1 = ctg4heapArray[nodeCounter].ctgID;

		if ( node1 == last_node )
		{
			node2 = last_cnt->contigID;
			bal_nd1 = getTwinCtg ( node1 );
			bal_nd2 = getTwinCtg ( node2 );
			temp_cnt = last_cnt;
			bal_cnt = checkConnect ( bal_nd2, bal_nd1 );
			next_cnt = checkConnect ( getTwinCtg ( ctg4heapArray[nodeCounter].ctgID ),
			                          getTwinCtg ( ctg4heapArray[nodeCounter - 1].ctgID ) );
			prev_cnt = checkConnect ( ctg4heapArray[nodeCounter - 1].ctgID, ctg4heapArray[nodeCounter].ctgID );
			setNextInScaf ( prev_cnt, temp_cnt );
			setNextInScaf ( bal_cnt, next_cnt );
			setPrevInScaf ( next_cnt, 1 );
		}
		else
		{
			next_node = last_cnt->contigID;
			bal_next_node = getTwinCtg ( next_node );
			j = nodeCounter;
			node1 = ctg4heapArray[j - 1].ctgID;
			node2 = ctg4heapArray[j].ctgID;
			ite_cnt = getCntBetween ( node2, next_node );
			bal_cnt = getCntBetween ( bal_next_node, getTwinCtg ( node2 ) );

			while ( !ite_cnt && node1 != last_node )
			{
				tmp_cnt = getCntBetween ( node1, node2 );
				bal_cnt = getCntBetween ( getTwinCtg ( node2 ), getTwinCtg ( node1 ) );
				setNextInScaf ( tmp_cnt, NULL );
				setPrevInScaf ( tmp_cnt, 0 );
				tmp_cnt->deleted = 1;
				setNextInScaf ( bal_cnt, NULL );
				setPrevInScaf ( bal_cnt, 0 );
				bal_cnt->deleted = 1;
				--j;
				node1 = ctg4heapArray[j - 1].ctgID;
				node2 = ctg4heapArray[j].ctgID;
				ite_cnt = getCntBetween ( node2, next_node );
				bal_cnt = getCntBetween ( bal_next_node, getTwinCtg ( node2 ) );
			}

			if ( !ite_cnt )
			{
				tmp_cnt = getCntBetween ( node1, node2 );
				gap = last_cnt->gapLen - tmp_cnt->gapLen - contig_array[node2].length;
				ite_cnt = allocateCN ( next_node, gap );
				ite_cnt->weight = 0;

				if ( cntLookupTable )
				{
					putCnt2LookupTable ( node2, ite_cnt );
				}

				ite_cnt->next = contig_array[node2].downwardConnect;
				contig_array[node2].downwardConnect = ite_cnt;
				bal_cnt = allocateCN ( getTwinCtg ( node2 ), gap );
				bal_cnt->weight = 0;

				if ( cntLookupTable )
				{
					putCnt2LookupTable ( bal_next_node, bal_cnt );
				}

				bal_cnt->next = contig_array[bal_next_node].downwardConnect;
				contig_array[bal_next_node].downwardConnect = bal_cnt;
			}

			ite_cnt->deleted = 0;
			ite_cnt->mask = 0;
			bal_cnt->deleted = 0;
			bal_cnt->mask = 0;
			tmp_cnt = getCntBetween ( node1, node2 );
			setNextInScaf ( tmp_cnt, ite_cnt );
			setPrevInScaf ( ite_cnt, 1 );
			tmp_cnt = getCntBetween ( getTwinCtg ( node2 ), getTwinCtg ( node1 ) );
			setNextInScaf ( bal_cnt, tmp_cnt );
			setPrevInScaf ( tmp_cnt, 1 );

			if ( last_cnt->nextInScaf )
			{
				setNextInScaf ( ite_cnt, last_cnt->nextInScaf );
				affected_node1 = last_cnt->nextInScaf->contigID;
				bal_affected_node1 = getTwinCtg ( affected_node1 );
				bal_affected_cnt = getCntBetween ( bal_affected_node1, bal_next_node );
				setNextInScaf ( bal_affected_cnt, bal_cnt );
				setPrevInScaf ( bal_cnt, 1 );
			}

			setNextInScaf ( last_cnt, NULL );
			setPrevInScaf ( last_cnt, 0 );
			last_cnt->deleted = 1;
			tmp_cnt = getCntBetween ( bal_next_node, bal_last_node );
			setNextInScaf ( tmp_cnt, NULL );
			setPrevInScaf ( tmp_cnt, 0 );
			tmp_cnt->deleted = 1;
		}
	}
}

/*************************************************
Function:
    checkOverlapInBetween_general
Description:
    Checks if adjacent contigs in the array have reasonable overlap.
Input:
    1. tolerance:       max percentage that overlap length accounts for
Output:
    None.
Return:
    1 if the overlap situation was resonable.
*************************************************/
boolean checkOverlapInBetween_general ( double tolerance )
{
	int i, gap;
	unsigned int node1, node2;
	int lenSum, lenOlp;
	CONNECT * cnt;
	lenSum = lenOlp = 0;

	for ( i = 1; i <= nodeCounter; i++ )
	{
		node1 = ctg4heapArray[i].ctgID;
		lenSum += contig_array[node1].length;
	}

	if ( lenSum < 1 )
		{ return 0; }

	for ( i = 1; i < nodeCounter; i++ )
	{
		node2 = ctg4heapArray[i + 1].ctgID;
		gap = ctg4heapArray[i + 1].dis - ctg4heapArray[i].dis
		      - contig_array[node2].length;

		if ( -gap > 0 )
		{
			node1 = ctg4heapArray[i].ctgID;
			cnt = getCntBetween ( node1, node2 );

			if ( cnt && cnt->gapLen > gap )
			{
				continue;
			}
			else if ( ( cnt = getCntBetween ( node2, node1 ) ) != NULL
			          && cnt->gapLen > gap )
			{
				continue;
			}
			else if ( -gap < overlaplen )
			{
				continue;
			}

			lenOlp += -gap;
		}

		if ( ( double ) lenOlp / lenSum > tolerance )
			{ return 0; }
	}

	return 1;
}

int canexchange = 0, exchange_num = 0, failexchange = 0;

/*************************************************
Function:
    checkConflictCnt_general
Description:
    Checks if the conflict connections between adjacent contigs
    in the array are acceptable.
Input:
    1. tolerance:       max percentage that conflict connections accounts for
Output:
    None.
Return:
    0 if acceptable.
*************************************************/
static boolean checkConflictCnt_general ( double tolerance )
{
	int i, j, gap;
	int supportCounter = 0;
	int objectCounter = 0;
	CONNECT * cnt;

	for ( i = 1; i < nodeCounter; i++ )
	{
		for ( j = i + 1; j <= nodeCounter; j++ )
		{
			cnt = checkConnect ( ctg4heapArray[i].ctgID, ctg4heapArray[j].ctgID );

			if ( cnt )
				{ supportCounter += cnt->weight; }

			cnt = checkConnect ( ctg4heapArray[j].ctgID, ctg4heapArray[i].ctgID );

			if ( cnt )
			{
				gap = ctg4heapArray[j].dis - ctg4heapArray[i].dis - contig_array[ctg4heapArray[j].ctgID].length;

				if ( gap > -overlaplen && gap >= cnt->gapLen && cnt->inherit == 0 )
				{
					objectCounter += cnt->weight;
				}
			}
		}
	}

	if ( supportCounter < 1 )
		{ return 1; }

	if ( ( double ) objectCounter / supportCounter < tolerance )
		{ return 0; }

	return 1;
}

/*************************************************
Function:
    getIndexInSortedSubgraph
Description:
    Gets contig's index in sorted array.
Input:
    1. node:        contig
    2. count:       array size
Output:
    None.
Return:
    Contig's index if found.
    -1 otherwise.
*************************************************/
static int getIndexInSortedSubgraph ( unsigned int node, int count )
{
	int index;

	for ( index = 0; index < count; ++index )
	{
		if ( nodesInSubInOrder[index] == node )
			{ return index; }
	}

	return -1;
}

/*************************************************
Function:
    getIndexInSortedSubgraph
Description:
    Gets contig's arc if contig has only one arc with weight > 1.
Input:
    1. node:        contig
Output:
    None.
Return:
    Pointer to arc if existed.
    NULL otherwise.
*************************************************/
static preARC * getValidArc ( unsigned int node )
{
	int num = 0;
	preARC * arc = contig_array[node].arcs;

	while ( arc )
	{
		if ( arc->multiplicity > 1 )
		{
			++num;

			if ( num > 1 )
			{
				return NULL;
			}
		}

		arc = arc->next;
	}

	return arc;
}

static boolean clearUpSubgraph()
{
	unsigned int i, ctg1, bal_ctg1, ctg2, bal_ctg2;
	int j, arc_num, num5, num3, index = 0, count = 0;
	preARC * arc;
	CONNECT * cnt;

	//put all contigs in "nodesInSub" array
	for ( i = 1; i <= nodeCounter; ++i )
	{
		nodesInSub[i - 1] = ctg4heapArray[i].ctgID;
	}

	for ( i = 0; i < nodeCounter; ++i )
	{
		ctg1 = nodesInSub[i];
		index = getIndexInSortedSubgraph ( ctg1, count );

		if ( index >= 0 && index < count - 1 ) //this contig is already in array
			{ continue; }

		bal_ctg1 = getTwinCtg ( ctg1 );
		num5 = 0;
		num3 = 0;
		arc_num = 0;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg1;
		arc = getValidArc ( ctg1 );

		while ( arc )
		{
			ctg2 = arc->to_ed;
			bal_ctg2 = getTwinCtg ( ctg2 );

			if ( ( arc = getValidArc ( bal_ctg2 ) ) == NULL )
			{
				break;
			}
			else if ( arc->to_ed != bal_ctg1 )
			{
				break;
			}

			ctg1 = ctg2;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg1;
			arc = getValidArc ( ctg1 );
		}

		ctg1 = nodesInSub[i];
		arc = getValidArc ( bal_ctg1 );

		while ( arc )
		{
			bal_ctg2 = arc->to_ed;
			ctg2 = getTwinCtg ( bal_ctg2 );

			if ( ( arc = getValidArc ( ctg2 ) ) == NULL )
			{
				break;
			}
			else if ( arc->to_ed != ctg1 )
			{
				break;
			}

			ctg1 = ctg2;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = ctg1;
			arc = getValidArc ( bal_ctg2 );
		}

		for ( j = num3 - 1; j >= 0; --j )
		{
			nodesInSubInOrder[index++] = * ( unsigned int * ) darrayGet ( scaf3, j );
		}

		for ( j = 0; j < num5; ++j )
		{
			nodesInSubInOrder[index++] = * ( unsigned int * ) darrayGet ( scaf5, j );
		}
	}
}

/*************************************************
Function:
    transferCnt2RemainNode
Description:
    In bubble structure, transfers connections of contig which is to
    be masked, to remained contig.
Input:
    1. maskNode:            contig to be masked
    2. remainNode:      remained contig
Output:
    None.
Return:
    None.
*************************************************/
static void transferCnt2RemainNode ( unsigned int maskNode, unsigned int remainNode )
{
	CONNECT * cnt = contig_array[maskNode].downwardConnect;
	CONNECT * bal_cnt, *nextCnt, *bal_nextCnt, *tmpCnt, *bal_tmpCnt;
	unsigned int nextNode1, bal_nextNode1, nextNode2, bal_nextNode2;
	unsigned int bal_maskNode = getTwinCtg ( maskNode ), bal_remainNode = getTwinCtg ( remainNode );
	int gap, weight, inherit;

	while ( cnt )
	{
		if ( cnt->mask )
		{
			cnt = cnt->next;
			continue;
		}

		nextNode1 = cnt->contigID;
		bal_nextNode1 = getTwinCtg ( nextNode1 );
		bal_cnt = getCntBetween ( bal_nextNode1, bal_maskNode );
		gap = cnt->gapLen;
		weight = cnt->weight;
		tmpCnt = getCntBetween ( remainNode, nextNode1 );

		if ( tmpCnt )
		{
			inherit = 0;
		}
		else
		{
			inherit = 1;
		}

		if ( cnt->nextInScaf )
		{
			nextNode2 = cnt->nextInScaf->contigID;
			bal_nextNode2 = getTwinCtg ( nextNode2 );
			nextCnt = getCntBetween ( nextNode1, nextNode2 );
			bal_nextCnt = getCntBetween ( bal_nextNode2, bal_nextNode1 );

			if ( nextNode1 != remainNode && nextNode2 != remainNode )
			{
				tmpCnt = add1Connect ( remainNode, nextNode1, gap, weight, inherit );
				bal_tmpCnt = add1Connect ( bal_nextNode1, bal_remainNode, gap, weight, inherit );
				tmpCnt->nextInScaf = nextCnt;
				tmpCnt->mask = 0;
				tmpCnt->deleted = 0;
				bal_nextCnt->nextInScaf = bal_tmpCnt;
				bal_tmpCnt->prevInScaf = 1;
				bal_tmpCnt->mask = 0;
				bal_tmpCnt->deleted = 0;
			}
			else
			{
				nextCnt->prevInScaf = 0;
				bal_nextCnt->nextInScaf = NULL;
			}
		}

		cnt->nextInScaf = NULL;
		cnt->prevInScaf = 0;
		cnt->mask = 1;
		cnt->deleted = 1;
		bal_cnt->nextInScaf = NULL;
		bal_cnt->prevInScaf = 0;
		bal_cnt->mask = 1;
		bal_cnt->deleted = 1;
		cnt = cnt->next;
	}
}

/*************************************************
Function:
    maskNodeCnt
Description:
    Masks contig's downstream connections in scaffold.
Input:
    1. node:        contig whose downstream connections are going to be masked
Output:
    None.
Return:
    None.
*************************************************/
static void maskNodeCnt ( unsigned int node )
{
	CONNECT * cnt = contig_array[node].downwardConnect;
	CONNECT * bal_cnt, *bal_nextCnt;
	unsigned int bal_nextNode1, bal_nextNode2;

	while ( cnt )
	{
		bal_nextNode1 = getTwinCtg ( cnt->contigID );
		bal_cnt = getCntBetween ( bal_nextNode1, getTwinCtg ( node ) );

		if ( cnt->nextInScaf )
		{
			cnt->nextInScaf->prevInScaf = 0;
			bal_nextNode2 = getTwinCtg ( cnt->nextInScaf->contigID );
			bal_nextCnt = getCntBetween ( bal_nextNode2, bal_nextNode1 );

			if ( !bal_nextCnt )
			{
				exit ( 1 );
			}

			bal_nextCnt->nextInScaf = NULL;
		}

		cnt->nextInScaf = NULL;
		cnt->prevInScaf = 0;
		cnt->mask = 1;
		cnt->deleted = 1;
		bal_cnt->nextInScaf = NULL;
		bal_cnt->prevInScaf = 0;
		bal_cnt->mask = 1;
		bal_cnt->deleted = 1;
		cnt = cnt->next;
	}
}

/*************************************************
Function:
    getEndKmers
Description:
    Gets contig's first and last kmers' sequences.
Input:
    1. seq:     contig sequence
    2. len:     contig length
    3. rev:     indicator of reversed complement
Output:
    1. firstKmer:       first kmer's sequence
    2. lastKmer:        last kmer's sequence
Return:
    None.
*************************************************/
static void getEndKmers ( char * seq, int len, int rev, char * firstKmer, char * lastKmer )
{
	int j;

	if ( 0 == rev )
	{
		for ( j = 0; j < overlaplen; ++j )
		{
			firstKmer[j] = int2base ( ( int ) getCharInTightString ( seq, j ) );
			lastKmer[j] = int2base ( ( int ) getCharInTightString ( seq, len - j - 1 ) );
		}
	}
	else
	{
		for ( j = 0; j < overlaplen; ++j )
		{
			firstKmer[j] = int2compbase ( ( int ) getCharInTightString ( seq, len - j - 1 ) );
			lastKmer[j] = int2compbase ( ( int ) getCharInTightString ( seq, j ) );
		}
	}
}

/*************************************************
Function:
    output_ctg
Description:
    Outputs contig sequence.
Input:
    1. ctg:     contig
    2. fo:      output file
Output:
    None.
Return:
    None.
*************************************************/
static void output_ctg ( unsigned int ctg, FILE * fo )
{
	if ( contig_array[ctg].length < 1 )
		{ return; }

	int len;
	unsigned int bal_ctg = getTwinCtg ( ctg );
	len = contig_array[ctg].length + overlaplen;
	int col = 0;

	if ( contig_array[ctg].seq )
	{
		fprintf ( fo, ">C%d %4.1f\n", ctg, ( double ) contig_array[ctg].cvg );
		outputTightStr ( fo, contig_array[ctg].seq, 0, len, len, 0, &col );
	}
	else if ( contig_array[bal_ctg].seq )
	{
		fprintf ( fo, ">C%d %4.1f\n", bal_ctg, ( double ) contig_array[ctg].cvg );
		outputTightStr ( fo, contig_array[bal_ctg].seq, 0, len, len, 0, &col );
	}

	fprintf ( fo, "\n" );
}


/*************************************************
Function:
    removeBubbleCtg
Description:
    Searchs bubble structure in sub-graph. Removes the contig
    with lower coverage and transfers its connection relation to
    the remained contig if found structure accords with following
    criterions:
    1. Both contigs' coverage smaller than a cutoff.
    2. The first and last kmers of two contigs are the same.
Input:
    None.
Output:
    None.
Return:
    Handled bubble structure number.
*************************************************/
static int removeBubbleCtg()
{
	int i, j, count, gap, SnpCounter = 0, conflict = 0;
	unsigned int node1, node2, bal_node1, bal_node2;
	int len1, len2, addLast = 0;
	char * tightStr1, *tightStr2;
	char firstKmer1[overlaplen + 1], lastKmer1[overlaplen + 1], firstKmer2[overlaplen + 1], lastKmer2[overlaplen + 1];
	CONNECT * cnt, *bal_cnt;
	count = 0;

	for ( i = 1; i < nodeCounter; ++i )
	{
		node1 = ctg4heapArray[i].ctgID;
		node2 = ctg4heapArray[i + 1].ctgID;
		bal_node1 = getTwinCtg ( node1 );
		bal_node2 = getTwinCtg ( node2 );
		cnt = getCntBetween ( node1, node2 );
		bal_cnt = getCntBetween ( node2, node1 );
		gap = ctg4heapArray[i + 1].dis - ctg4heapArray[i].dis - ( int ) contig_array[node2].length;

		if ( gap >= 0 || contig_array[node1].cvg >= cvg4SNP || contig_array[node2].cvg >= cvg4SNP || cnt || bal_cnt )
		{
			nodesInSubInOrder[count] = node1;
			nodeDistanceInOrder[count++] = ctg4heapArray[i].dis;
			continue;
		}

		len1 = contig_array[node1].length + overlaplen;
		len2 = contig_array[node2].length + overlaplen;

		if ( contig_array[node1].seq )
		{
			getEndKmers ( contig_array[node1].seq, len1, 0, firstKmer1, lastKmer1 );
		}
		else
		{
			getEndKmers ( contig_array[bal_node1].seq, len1, 1, firstKmer1, lastKmer1 );
		}

		if ( contig_array[node2].seq )
		{
			getEndKmers ( contig_array[node2].seq, len2, 0, firstKmer2, lastKmer2 );
		}
		else
		{
			getEndKmers ( contig_array[bal_node2].seq, len2, 1, firstKmer2, lastKmer2 );
		}

		for ( j = 0; j < overlaplen; ++j )
		{
			if ( firstKmer1[j] != firstKmer2[j] || lastKmer1[j] != lastKmer2[j] )
			{
				nodesInSubInOrder[count] = node1;
				nodeDistanceInOrder[count++] = ctg4heapArray[i].dis;
				conflict = 1;
				break;
			}
		}

		if ( 1 == conflict )
		{
			conflict = 0;
			continue;
		}

		++SnpCounter;

		if ( contig_array[node1].bubbleInScaff == 0 || contig_array[node2].bubbleInScaff == 0 )
		{
			contig_array[node1].bubbleInScaff = 1;
			contig_array[bal_node1].bubbleInScaff = 1;
			contig_array[node2].bubbleInScaff = 1;
			contig_array[bal_node2].bubbleInScaff = 1;
			output_ctg ( node1, snp_fp );
			output_ctg ( node2, snp_fp );
		}

		if ( contig_array[node1].cvg > contig_array[node2].cvg || ( len1 > len2 && contig_array[node1].cvg == contig_array[node2].cvg ) )
		{
			if ( i == nodeCounter - 1 )
			{
				nodesInSubInOrder[count] = node1;
				nodeDistanceInOrder[count++] = ctg4heapArray[i].dis;
				addLast = 1;
			}

			transferCnt2RemainNode ( node2, node1 );
			transferCnt2RemainNode ( bal_node2, bal_node1 );
			contig_array[node2].mask = 1;
			contig_array[bal_node2].mask = 1;
			ctg4heapArray[i + 1].ctgID = node1;
			ctg4heapArray[i + 1].dis = ctg4heapArray[i].dis;
		}
		else
		{
			if ( i == nodeCounter - 1 )
			{
				nodesInSubInOrder[count] = node2;
				nodeDistanceInOrder[count++] = ctg4heapArray[i + 1].dis;
				addLast = 1;
			}

			transferCnt2RemainNode ( node1, node2 );
			transferCnt2RemainNode ( bal_node1, bal_node2 );
			contig_array[node1].mask = 1;
			contig_array[getTwinCtg ( node1 )].mask = 1;
		}
	}

	if ( 0 == addLast )
	{
		nodesInSubInOrder[count] = ctg4heapArray[nodeCounter].ctgID;
		nodeDistanceInOrder[count++] = ctg4heapArray[nodeCounter].dis;
	}

	for ( i = 0; i < count; ++i )
	{
		ctg4heapArray[i + 1].ctgID = nodesInSubInOrder[i];
		ctg4heapArray[i + 1].dis = nodeDistanceInOrder[i];
	}

	nodeCounter = count;
	return SnpCounter;
}

/*************************************************
Function:
    general_linearization
Description:
    Picks up sub-graph and trys to line involved contigs.
Input:
    1. strict:      indicator of whether using strict cutoff to deal with sub-graph
Output:
    None.
Return:
    None.
*************************************************/
static void general_linearization ( boolean strict )
{
	unsigned int i;
	int subCounter = 0;
	int out_num;
	boolean flag;
	int conflCounter = 0, overlapCounter = 0, eligibleCounter = 0;
	int SNPCtgCounter = 0;
	double overlapTolerance, conflTolerance;
	canexchange = 0, exchange_num = 0, failexchange = 0;
	fprintf ( stderr, "Start to linearize sub-graph.\n" );

	for ( i = num_ctg; i > 0; i-- )
	{
		if ( contig_array[i].mask )
			{ continue; }

		out_num = validConnect ( i, NULL );

		if ( out_num < 2 )
			{ continue; }

		flag = pickUpGeneralSubgraph ( i, MaxNodeInSub );

		if ( !flag )
		{
			continue;
		}

		subCounter++;
		qsort ( &ctg4heapArray[1], nodeCounter, sizeof ( CTGinHEAP ), cmp_ctg );

		if ( Insert_size < 1000 && cvg4SNP > 0.001 )
		{
			SNPCtgCounter += removeBubbleCtg();
		}

		flag = checkEligible();

		if ( !flag )
		{
			eligibleCounter++;
			setInGraph ( 0 );
			continue;
		}

		if ( strict )
		{
			overlapTolerance = OverlapPercent;
			conflTolerance = ConflPercent;
		}
		else
		{
			overlapTolerance = 2 * OverlapPercent;
			conflTolerance = 2 * ConflPercent;
		}

		flag = checkOverlapInBetween_general ( overlapTolerance );

		if ( !flag )
		{
			overlapCounter++;
			setInGraph ( 0 );
			continue;
		}

		flag = checkConflictCnt_general ( conflTolerance );

		if ( flag )
		{
			conflCounter++;
			setInGraph ( 0 );
			continue;
		}

		arrangeNodes_general();
		setInGraph ( 0 );
	}

	fprintf ( stderr, " Picked sub-graphs          %d\n Connection-conflict        %d\n Significant overlapping    %d\n Eligible                   %d\n Bubble structures          %d\n", subCounter, conflCounter, overlapCounter, eligibleCounter, SNPCtgCounter );
}

/****       the fowllowing codes for detecting and break down scaffold at weak point  **********/

/*************************************************
Function:
    smallScaf
Description:
    Counts constructed scaffold number by using short insert size
    paired-end reads and sets connections' flags: "bySmall" and
    "smallIns".
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void smallScaf()
{
	unsigned int i, ctg, bal_ctg, prevCtg;
	int counter = 0;
	CONNECT * bindCnt, *cnt;

	for ( i = 1; i <= num_ctg; i++ )
		{ contig_array[i].flag = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].flag || contig_array[i].mask || !contig_array[i].downwardConnect )
			{ continue; }

		bindCnt = getBindCnt ( i );

		if ( !bindCnt )
			{ continue; }

		counter++;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg ( i )].flag = 1;
		prevCtg = getTwinCtg ( i );

		while ( bindCnt )
		{
			ctg = bindCnt->contigID;
			bal_ctg = getTwinCtg ( ctg );
			bindCnt->bySmall = 1;
			cnt = getCntBetween ( bal_ctg, prevCtg );

			if ( cnt )
				{ cnt->bySmall = 1; }

			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCtg = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg ( i );
		bindCnt = getBindCnt ( ctg );
		prevCtg = i;

		while ( bindCnt )
		{
			ctg = bindCnt->contigID;
			bal_ctg = getTwinCtg ( ctg );
			bindCnt->bySmall = 1;
			cnt = getCntBetween ( bal_ctg, prevCtg );

			if ( cnt )
				{ cnt->bySmall = 1; }

			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCtg = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}
	}

	fprintf ( stderr, "Report from smallScaf: %d scaffolds by smallPE.\n", counter );

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( !contig_array[i].downwardConnect )
			{ continue; }

		cnt = contig_array[i].downwardConnect;

		while ( cnt )
		{
			cnt->smallIns = 1;
			cnt = cnt->next;
		}
	}
}

/*************************************************
Function:
    clearNewInsFlag
Description:
    Clears conections' flag "newIns".
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void clearNewInsFlag()
{
	int i = 1;
	CONNECT * cnt;

	for ( ; i <= num_ctg; i++ )
	{
		cnt = contig_array[i].downwardConnect;

		while ( cnt )
		{
			cnt->newIns = 0;

			if ( Insert_size > 15000 )
			{
				cnt = cnt->next;
				continue;
			}

			cnt->maxGap = 0;
			cnt = cnt->next;
		}
	}
}

/*************************************************
Function:
    putItem2Sarray
Description:
    Checks whether a scaffold ID is in array. If not, puts it in to array.
Input:
    1. scaf:            scaffold ID to be checked and put in to array
    2. wt:          weight of connection to scaffold ID
    3. SCAF:            array of scaffold IDs
    4. WT:          array of connections' weight
    5. counter:     size of scaffold ID array
Output:
    1. SCAF:        updated array of scaffold IDs
    2. WT:      updated array of connections' weight
Return:
    1 if operation of putting succeeded.
*************************************************/
static boolean putItem2Sarray ( unsigned int scaf, int wt, DARRAY * SCAF, DARRAY * WT, int counter )
{
	int i;
	unsigned int * scafP, *wtP;

	for ( i = 0; i < counter; i++ )
	{
		scafP = ( unsigned int * ) darrayGet ( SCAF, i );

		if ( ( *scafP ) == scaf )
		{
			wtP = ( unsigned int * ) darrayGet ( WT, i );
			*wtP = ( *wtP + wt );
			return 0;
		}
	}

	scafP = ( unsigned int * ) darrayPut ( SCAF, counter );
	wtP = ( unsigned int * ) darrayPut ( WT, counter );
	*scafP = scaf;
	*wtP = wt;
	return 1;
}

/*************************************************
Function:
    getDSLink2Scaf
Description:
    Gets downstream connections to other scaffolds and stores
    related information including other scaffolds' ID and connections'
    weight.
Input:
    1. scafStack:       scaffold stack consisting of contigs
    2. SCAF:            array to store other scaffolds' IDs
    3. WT:          array to store connections' weight
    4. total_len:       length sum of contigs used to search connections
                    to other scaffolds
Output:
    1. SCAF:        updated array of scaffold IDs
    2. WT:      updated array of connections' weight
Return:
    Number of other scaffolds being connected.
*************************************************/
static int getDSLink2Scaf ( STACK * scafStack, DARRAY * SCAF, DARRAY * WT, int total_len )
{
	CONNECT * ite_cnt;
	CONNECT * bind_cnt;
	unsigned int ctg, targetCtg, bal_targetCtg, *pt;
	int counter = 0;
	int len = 0, gap;
	boolean inc;
	stackRecover ( scafStack );

	while ( ( pt = ( unsigned int * ) stackPop ( scafStack ) ) != NULL )
	{
		ctg = *pt;
		bind_cnt = getBindCnt ( ctg );
		gap = bind_cnt ? bind_cnt->gapLen : 0;
		len += contig_array[ctg].length + gap;

		if ( ( contig_array[ctg].mask && contig_array[ctg].length < 500 ) || !contig_array[ctg].downwardConnect
		        || total_len - len > Insert_size )
		{
			continue;
		}

		ite_cnt = contig_array[ctg].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->newIns != 1 )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			targetCtg = ite_cnt->contigID;
			bal_targetCtg = getTwinCtg ( targetCtg );

			if ( ( ite_cnt->mask && contig_array[targetCtg].length < 500 ) || ite_cnt->singleInScaf
			        || ite_cnt->nextInScaf || ite_cnt->prevInScaf || ite_cnt->inherit )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			if ( contig_array[ctg].from_vt == contig_array[targetCtg].from_vt // on the same scaff
			        || ( targetCtg == contig_array[targetCtg].from_vt && bal_targetCtg == contig_array[bal_targetCtg].from_vt ) ) //targetCtg isn't in any scaffold
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			inc = putItem2Sarray ( contig_array[targetCtg].from_vt, ite_cnt->weight, SCAF, WT, counter );

			if ( inc )
				{ counter++; }

			ite_cnt = ite_cnt->next;
		}
	}

	return counter;
}

/*************************************************
Function:
    getScaffold
Description:
    Starting from a contig, gets downstream contig chain in
    a scaffold.
Input:
    1. start:       start contig
Output:
    1. scafStack:       stack to store contig chain
Return:
    Length of contig chain.
*************************************************/
static int getScaffold ( unsigned int start, STACK * scafStack )
{
	int len = contig_array[start].length;
	unsigned int * pt, ctg;
	emptyStack ( scafStack );
	pt = ( unsigned int * ) stackPush ( scafStack );
	*pt = start;
	CONNECT * bindCnt = getBindCnt ( start );

	while ( bindCnt )
	{
		ctg = bindCnt->contigID;
		pt = ( unsigned int * ) stackPush ( scafStack );
		*pt = ctg;
		len += bindCnt->gapLen + contig_array[ctg].length;
		bindCnt = bindCnt->nextInScaf;
	}

	stackBackup ( scafStack );
	return len;
}

/*************************************************
Function:
    isLinkReliable
Description:
    Checks whether there is a reliable connection.
Input:
    1. WT:      weight of connections
    2. count:       connection number
Output:
    None.
Return:
    1 if a reliable connection was found.
*************************************************/
static boolean isLinkReliable ( DARRAY * WT, int count )
{
	int i;

	for ( i = 0; i < count; i++ )
		if ( * ( int * ) darrayGet ( WT, i ) >= weakPE )
			{ return 1; }

	return 0;
}

/*************************************************
Function:
    getWtFromSarray
Description:
    Gets weight of connection to specified scaffold.
Input:
    1. SCAF:        array of scaffold ID
    2. WT:      array of onnection weight
    3. count:       array size
    4. scaf:        specified scaffold
Output:
    None.
Return:
    Weight of connection was found.
    0 otherwise.
*************************************************/
static int getWtFromSarray ( DARRAY * SCAF, DARRAY * WT, int count, unsigned int scaf )
{
	int i;

	for ( i = 0; i < count; i++ )
		if ( * ( unsigned int * ) darrayGet ( SCAF, i ) == scaf )
			{ return * ( int * ) darrayGet ( WT, i ); }

	return 0;
}

/*************************************************
Function:
    switch2twin
Description:
    Transfers contigs to their reversed complements.
Input:
    1. scafStack:       stack of contigs
Output:
    1. scafStack:       updated stack of contigs
Return:
    None.
*************************************************/
static void switch2twin ( STACK * scafStack )
{
	unsigned int * pt;
	stackRecover ( scafStack );

	while ( ( pt = ( unsigned int * ) stackPop ( scafStack ) ) != NULL )
		{ *pt = getTwinCtg ( *pt ); }
}

/*************************************************
Function:
    recoverLinks
Description:
    Recovers contigs' connections to other scaffold if they
    accords with some criterions.
Input:
    1. scafStack:       stack of contigs
Output:
    None.
Return:
    None.
*************************************************/
static void recoverLinks ( STACK * scafStack )
{
	CONNECT * ite_cnt;
	unsigned int ctg, targetCtg, *pt;
	int counter = 0;
	boolean inc;
	unsigned int bal_ctg;
	stackRecover ( scafStack );

	while ( ( pt = ( unsigned int * ) stackPop ( scafStack ) ) != NULL )
	{
		ctg = *pt;

		if ( contig_array[ctg].mask || !contig_array[ctg].downwardConnect )
		{
			continue;
		}

		ite_cnt = contig_array[ctg].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->mask || ite_cnt->singleInScaf || ite_cnt->nextInScaf || ite_cnt->prevInScaf || ite_cnt->inherit || ite_cnt->weight < weakPE )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			targetCtg = ite_cnt->contigID;

			if ( contig_array[ctg].from_vt == contig_array[targetCtg].from_vt )     // on the same scaff
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			setConnectDelete ( ctg, targetCtg, 0, 0 );
			ite_cnt = ite_cnt->next;
		}

		bal_ctg = getTwinCtg ( ctg );

		if ( contig_array[bal_ctg].mask || !contig_array[bal_ctg].downwardConnect )
		{
			continue;
		}

		ite_cnt = contig_array[bal_ctg].downwardConnect;

		while ( ite_cnt )
		{
			if ( ite_cnt->mask || ite_cnt->singleInScaf || ite_cnt->nextInScaf || ite_cnt->prevInScaf || ite_cnt->inherit || ite_cnt->weight < weakPE )
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			targetCtg = ite_cnt->contigID;

			if ( contig_array[bal_ctg].from_vt == contig_array[targetCtg].from_vt )     // on the same scaff
			{
				ite_cnt = ite_cnt->next;
				continue;
			}

			setConnectDelete ( bal_ctg, targetCtg, 0, 0 );
			ite_cnt = ite_cnt->next;
		}
	}
}
/*
           ------>
 scaf1 --- --- -- ---
                         scaf2 -- --- --- --
                                   ---->
*/
/*************************************************
Function:
    checkScafConsist
Description:
    Checks  if both sets of contigs in two stacks have reliable
    connection to different scaffolds.
Input:
    1. scafStack1:      stack of contig set one
    2. len1:            total length of contigs in set one
    3.  scafStack2: stack of contig set two
    4. len2:            total length of contigs in set two
Output:
    None.
Return:
    0 if both sets of contigs did not have reliable connection
    to different scaffolds.
*************************************************/
static boolean checkScafConsist ( STACK * scafStack1, int len1, STACK * scafStack2, int len2 )
{
	DARRAY * downwardTo1 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) ); // scaf links to those scaffolds
	DARRAY * downwardTo2 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	DARRAY * downwardWt1 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) ); // scaf links to scaffolds with those wt
	DARRAY * downwardWt2 = ( DARRAY * ) createDarray ( 1000, sizeof ( unsigned int ) );
	int linkCount1 = getDSLink2Scaf ( scafStack1, downwardTo1, downwardWt1, len1 );
	int linkCount2 = getDSLink2Scaf ( scafStack2, downwardTo2, downwardWt2, len2 );

	if ( !linkCount1 )
	{
		freeDarray ( downwardTo1 );
		freeDarray ( downwardTo2 );
		freeDarray ( downwardWt1 );
		freeDarray ( downwardWt2 );
		return 1;
	}

	boolean flag1 = isLinkReliable ( downwardWt1, linkCount1 );

	if ( !flag1 )
	{
		freeDarray ( downwardTo1 );
		freeDarray ( downwardTo2 );
		freeDarray ( downwardWt1 );
		freeDarray ( downwardWt2 );
		return 1;
	}

	unsigned int scaf;
	int i, wt1, wt2, ret = 1;

	for ( i = 0; i < linkCount1; i++ )
	{
		wt1 = * ( int * ) darrayGet ( downwardWt1, i );

		if ( wt1 < weakPE )
			{ continue; }

		scaf = * ( unsigned int * ) darrayGet ( downwardTo1, i );
		wt2 = getWtFromSarray ( downwardTo2, downwardWt2, linkCount2, scaf );

		if ( wt2 < 1 )
		{
			ret = 0;
			break;
		}
	}

	if ( ret == 0 )
	{
		if ( linkCount1 && flag1 )
		{
			recoverLinks ( scafStack1 );
		}

		if ( linkCount2 )
		{
			recoverLinks ( scafStack2 );
		}
	}

	freeDarray ( downwardTo1 );
	freeDarray ( downwardTo2 );
	freeDarray ( downwardWt1 );
	freeDarray ( downwardWt2 );
	return ret;
}

/*************************************************
Function:
    setBreakPoints
Description:
    Sets break points of scaffold at weak connections.
Input:
    1. ctgArray:        array of contig in scaffold
    2. count:           total contig number
    3. weakest:     contig whose downstream connection was
                    estimated break point
Output:
    1. start:       the first contig with weak downstream connection
    2. finish:      the last contig with downstream connection
Return:
    None.
*************************************************/
static void setBreakPoints ( DARRAY * ctgArray, int count, int weakest,
                             int * start, int * finish )
{
	int index = weakest - 1;
	unsigned int thisCtg;
	unsigned int nextCtg = * ( unsigned int * ) darrayGet ( ctgArray, weakest );
	CONNECT * cnt;
	*start = weakest;

	while ( index >= 0 )
	{
		thisCtg = * ( unsigned int * ) darrayGet ( ctgArray, index );
		cnt = getCntBetween ( thisCtg, nextCtg );

		if ( cnt->maxGap > 2 )
			{ break; }
		else
			{ *start = index; }

		nextCtg = thisCtg;
		index--;
	}

	unsigned int prevCtg = * ( unsigned int * ) darrayGet ( ctgArray, weakest + 1 );
	*finish = weakest + 1;
	index = weakest + 2;

	while ( index < count )
	{
		thisCtg = * ( unsigned int * ) darrayGet ( ctgArray, index );
		cnt = getCntBetween ( prevCtg, thisCtg );

		if ( cnt->maxGap > 2 )
			{ break; }
		else
			{ *finish = index; }

		prevCtg = thisCtg;
		index++;
	}
}

/*************************************************
Function:
    changeScafEnd
Description:
    Changes "to_vt" of contigs in a scaffold to a specified contig.
Input:
    1. scafStack:       stack of contigs in scaffold
    2. end:         specified contig
Output:
    None.
Return:
    None.
*************************************************/
static void changeScafEnd ( STACK * scafStack, unsigned int end )
{
	unsigned int ctg, *pt;
	unsigned int start = getTwinCtg ( end );
	stackRecover ( scafStack );

	while ( ( pt = ( unsigned int * ) stackPop ( scafStack ) ) != NULL )
	{
		ctg = *pt;
		contig_array[ctg].to_vt = end;
		contig_array[getTwinCtg ( ctg )].from_vt = start;
	}
}

/*************************************************
Function:
    changeScafBegin
Description:
    Changes "from_vt" of contigs in a scaffold to a specified
    contig.
Input:
    1. scafStack:       stack of contigs in scaffold
    2. start:           specified contig
Output:
    None.
Return:
    None.
*************************************************/
static void changeScafBegin ( STACK * scafStack, unsigned int start )
{
	unsigned int ctg, *pt;
	unsigned int end = getTwinCtg ( start );
	stackRecover ( scafStack );

	while ( ( pt = ( unsigned int * ) stackPop ( scafStack ) ) != NULL )
	{
		ctg = *pt;
		contig_array[ctg].from_vt = start;
		contig_array[getTwinCtg ( ctg )].to_vt = end;
	}
}

/*************************************************
Function:
    detectBreakScaff
Description:
    Detects and breaks the weak connection between contigs
    in scaffold according to the support of large insert size
    paired-end reads.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void detectBreakScaf()
{
	unsigned int i, avgPE, scafLen, len, ctg, bal_ctg, prevCtg, thisCtg;
	long long peCounter, linkCounter;
	int num3, num5, weakPoint, tempCounter, j, t, counter = 0;
	CONNECT * bindCnt, *cnt, *weakCnt;
	STACK * scafStack1 = ( STACK * ) createStack ( 1000, sizeof ( unsigned int ) );
	STACK * scafStack2 = ( STACK * ) createStack ( 1000, sizeof ( unsigned int ) );

	for ( i = 1; i <= num_ctg; i++ )
		{ contig_array[i].flag = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].flag || contig_array[i].mask || !contig_array[i].downwardConnect )
			{ continue; }

		bindCnt = getBindCnt ( i );

		if ( !bindCnt )
			{ continue; }

		//first scan to get the average coverage by longer pe
		num5 = num3 = peCounter = linkCounter = 0;
		scafLen = contig_array[i].length;
		ctg = i;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg ( i )].flag = 1;

		while ( bindCnt )
		{
			if ( !bindCnt->bySmall )
				{ break; }

			linkCounter++;
			peCounter += bindCnt->maxGap;
			ctg = bindCnt->contigID;
			scafLen += contig_array[ctg].length;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg ( i );
		bindCnt = getBindCnt ( ctg );

		while ( bindCnt )
		{
			if ( !bindCnt->bySmall )
				{ break; }

			linkCounter++;
			peCounter += bindCnt->maxGap;
			ctg = bindCnt->contigID;
			scafLen += contig_array[ctg].length;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}

		if ( linkCounter < 1 || scafLen < 5000 )
			{ continue; }

		avgPE = peCounter / linkCounter;

		if ( avgPE < 10 )
			{ continue; }

		tempCounter = 0;

		for ( j = num3 - 1; j >= 0; j-- )
			* ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf3, j );

		for ( j = 0; j < num5; j++ )
			* ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf5, j );

		prevCtg = * ( unsigned int * ) darrayGet ( tempArray, 0 );
		weakCnt = NULL;
		weakPoint = 0;
		len = contig_array[prevCtg].length;

		for ( t = 1; t < tempCounter; t++ )
		{
			thisCtg = * ( unsigned int * ) darrayGet ( tempArray, t );

			if ( len < 2000 )
			{
				len += contig_array[thisCtg].length;
				prevCtg = thisCtg;
				continue;
			}
			else if ( len > scafLen - 2000 )
				{ break; }

			len += contig_array[thisCtg].length;

			if ( contig_array[prevCtg].from_vt != contig_array[thisCtg].from_vt ||
			        contig_array[prevCtg].indexInScaf > contig_array[thisCtg].indexInScaf )
			{
				prevCtg = thisCtg;
				continue;
			}

			cnt = getCntBetween ( prevCtg, thisCtg );

			if ( !weakCnt || weakCnt->maxGap > cnt->maxGap )
			{
				weakCnt = cnt;
				weakPoint = t;
			}

			prevCtg = thisCtg;
		}

		if ( !weakCnt || ( weakCnt->maxGap > 2 && weakCnt->maxGap > avgPE / 5 ) )
			{ continue; }

		prevCtg = * ( unsigned int * ) darrayGet ( tempArray, weakPoint - 1 );
		thisCtg = * ( unsigned int * ) darrayGet ( tempArray, weakPoint );

		if ( contig_array[prevCtg].from_vt != contig_array[thisCtg].from_vt ||
		        contig_array[prevCtg].indexInScaf > contig_array[thisCtg].indexInScaf )
		{
			fprintf ( stderr, "contig %d and %d not on the same scaff\n", prevCtg, thisCtg );
			continue;
		}

		setConnectWP ( prevCtg, thisCtg, 1 );
		int index1, index2;
		setBreakPoints ( tempArray, tempCounter, weakPoint - 1, &index1, &index2 );
		unsigned int start = * ( unsigned int * ) darrayGet ( tempArray, index1 );
		unsigned int finish = * ( unsigned int * ) darrayGet ( tempArray, index2 );
		int len1 = getScaffold ( getTwinCtg ( start ), scafStack1 );
		int len2 = getScaffold ( finish, scafStack2 );

		if ( len1 < 2000 || len2 < 2000 )
			{ continue; }

		switch2twin ( scafStack1 );
		int flag1 = checkScafConsist ( scafStack1, len1, scafStack2, len2 );
		switch2twin ( scafStack1 );
		switch2twin ( scafStack2 );
		int flag2 = checkScafConsist ( scafStack2, len2, scafStack1, len1 );

		if ( !flag1 || !flag2 )
		{
			changeScafBegin ( scafStack1, getTwinCtg ( start ) );
			changeScafEnd ( scafStack2, getTwinCtg ( finish ) );
			//unbind links
			unsigned int nextCtg = * ( unsigned int * ) darrayGet ( tempArray, index1 + 1 );
			thisCtg = * ( unsigned int * ) darrayGet ( tempArray, index1 );
			cnt = getCntBetween ( getTwinCtg ( nextCtg ), getTwinCtg ( thisCtg ) );

			if ( cnt->nextInScaf )
			{
				prevCtg = getTwinCtg ( cnt->nextInScaf->contigID );
				cnt->nextInScaf->prevInScaf = 0;
				cnt = getCntBetween ( prevCtg, thisCtg );
				cnt->nextInScaf = NULL;
			}

			prevCtg = * ( unsigned int * ) darrayGet ( tempArray, index2 - 1 );
			thisCtg = * ( unsigned int * ) darrayGet ( tempArray, index2 );
			cnt = getCntBetween ( prevCtg, thisCtg );

			if ( cnt->nextInScaf )
			{
				nextCtg = cnt->nextInScaf->contigID;
				cnt->nextInScaf->prevInScaf = 0;
				cnt = getCntBetween ( getTwinCtg ( nextCtg ), getTwinCtg ( thisCtg ) );
				cnt->nextInScaf = NULL;
			}

			prevCtg = * ( unsigned int * ) darrayGet ( tempArray, index1 );

			for ( t = index1 + 1; t <= index2; t++ )
			{
				thisCtg = * ( unsigned int * ) darrayGet ( tempArray, t );
				cnt = getCntBetween ( prevCtg, thisCtg );
				cnt->mask = 1;
				cnt->nextInScaf = NULL;
				cnt->prevInScaf = 0;
				cnt = getCntBetween ( getTwinCtg ( thisCtg ), getTwinCtg ( prevCtg ) );
				cnt->mask = 1;
				cnt->nextInScaf = NULL;
				cnt->prevInScaf = 0;
				prevCtg = thisCtg;
			}

			counter++;
		}
	}

	freeStack ( scafStack1 );
	freeStack ( scafStack2 );
	fprintf ( stderr, "Report from checkScaf: %d scaffold segments broken.\n", counter );
}


/*************************************************
Function:
    detectBreakScaff
Description:
    Detects and breaks the weak connection between contigs
    in scaffold according to the support of large insert size
    paired-end reads.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void detectBreakScaff()
{
	unsigned int i, avgPE, scafLen, len, ctg, bal_ctg, prevCtg, thisCtg;
	long long peCounter, linkCounter;
	int num3, num5, weakPoint, tempCounter, j, t, counter = 0;
	int newInsNum;
	CONNECT * bindCnt, *cnt, *weakCnt;
	CONNECT * bal_cnt;
	STACK * scafStack1 = ( STACK * ) createStack ( 1000, sizeof ( unsigned int ) );
	STACK * scafStack2 = ( STACK * ) createStack ( 1000, sizeof ( unsigned int ) );

	for ( i = 1; i <= num_ctg; i++ )
		{ contig_array[i].flag = 0; }

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].flag || contig_array[i].mask || !contig_array[i].downwardConnect )
			{ continue; }

		bindCnt = getBindCnt ( i );

		if ( !bindCnt )
			{ continue; }

		//first scan to get the average coverage by longer pe
		num5 = num3 = peCounter = linkCounter = 0;
		scafLen = contig_array[i].length;
		ctg = i;
		* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg ( i )].flag = 1;

		while ( bindCnt )
		{
			linkCounter++;
			peCounter += bindCnt->maxGap;
			ctg = bindCnt->contigID;
			scafLen += bindCnt->gapLen + contig_array[ctg].length;
			* ( unsigned int * ) darrayPut ( scaf5, num5++ ) = ctg;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg ( i );
		bindCnt = getBindCnt ( ctg );

		while ( bindCnt )
		{
			linkCounter++;
			peCounter += bindCnt->maxGap;
			ctg = bindCnt->contigID;
			scafLen += bindCnt->gapLen + contig_array[ctg].length;
			bal_ctg = getTwinCtg ( ctg );
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			* ( unsigned int * ) darrayPut ( scaf3, num3++ ) = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}

		if ( scafLen < Insert_size )
			{ continue; }

		avgPE = peCounter / linkCounter;

		if ( avgPE < 10 )
			{ continue; }

		tempCounter = 0;

		for ( j = num3 - 1; j >= 0; j-- )
			* ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf3, j );

		for ( j = 0; j < num5; j++ )
			* ( unsigned int * ) darrayPut ( tempArray, tempCounter++ ) =
			    * ( unsigned int * ) darrayGet ( scaf5, j );

		prevCtg = * ( unsigned int * ) darrayGet ( tempArray, 0 );
		weakCnt = NULL;
		weakPoint = 0;
		len = contig_array[prevCtg].length;

		for ( t = 1; t < tempCounter; t++ )
		{
			newInsNum = 0;
			thisCtg = * ( unsigned int * ) darrayGet ( tempArray, t );
			cnt = contig_array[thisCtg].downwardConnect;

			while ( cnt )
			{
				if ( cnt->newIns == 1 )
				{
					ctg = cnt->contigID;
					bal_ctg = getTwinCtg ( ctg );
					newInsNum++;
				}

				cnt = cnt->next;
			}

			bal_cnt = contig_array[getTwinCtg ( thisCtg )].downwardConnect;
			cnt = getCntBetween ( prevCtg, thisCtg );

			if ( len < Insert_size )
			{
				len += cnt->gapLen + contig_array[thisCtg].length;
				prevCtg = thisCtg;
				continue;
			}
			else if ( len > scafLen - Insert_size )
				{ break; }

			len += cnt->gapLen + contig_array[thisCtg].length;

			if ( contig_array[prevCtg].from_vt != contig_array[thisCtg].from_vt ||
			        contig_array[prevCtg].indexInScaf > contig_array[thisCtg].indexInScaf )
			{
				prevCtg = thisCtg;
				continue;
			}

			if ( !weakCnt || weakCnt->maxGap > cnt->maxGap )
			{
				weakCnt = cnt;
				weakPoint = t;
			}

			prevCtg = thisCtg;
		}

		if ( !weakCnt || ( weakCnt->maxGap > 2 && weakCnt->maxGap > avgPE / 5 ) )
			{ continue; }

		prevCtg = * ( unsigned int * ) darrayGet ( tempArray, weakPoint - 1 );
		thisCtg = * ( unsigned int * ) darrayGet ( tempArray, weakPoint );

		if ( contig_array[prevCtg].from_vt != contig_array[thisCtg].from_vt ||
		        contig_array[prevCtg].indexInScaf > contig_array[thisCtg].indexInScaf )
		{
			printf ( "contig %d and %d not on the same scaff\n", prevCtg, thisCtg );
			continue;
		}

		setConnectWP ( prevCtg, thisCtg, 1 );
		// set start and end to break down the scaffold
		int index1, index2;
		setBreakPoints ( tempArray, tempCounter, weakPoint - 1, &index1, &index2 );
		unsigned int start = * ( unsigned int * ) darrayGet ( tempArray, index1 );
		unsigned int finish = * ( unsigned int * ) darrayGet ( tempArray, index2 );
		int len1 = getScaffold ( getTwinCtg ( start ), scafStack1 );
		int len2 = getScaffold ( finish, scafStack2 );

		if ( len1 < Insert_size || len2 < Insert_size )
			{ continue; }

		switch2twin ( scafStack1 );
		int flag1 = checkScafConsist ( scafStack1, len1, scafStack2, len2 );
		switch2twin ( scafStack1 );
		switch2twin ( scafStack2 );
		int flag2 = checkScafConsist ( scafStack2, len2, scafStack1, len1 );

		if ( !flag1 || !flag2 )
		{
			changeScafBegin ( scafStack1, getTwinCtg ( start ) );
			changeScafEnd ( scafStack2, getTwinCtg ( finish ) );
			//unbind links
			unsigned int nextCtg = * ( unsigned int * ) darrayGet ( tempArray, index1 + 1 );
			thisCtg = * ( unsigned int * ) darrayGet ( tempArray, index1 );
			cnt = getCntBetween ( getTwinCtg ( nextCtg ), getTwinCtg ( thisCtg ) );

			if ( cnt->nextInScaf )
			{
				prevCtg = getTwinCtg ( cnt->nextInScaf->contigID );
				cnt->nextInScaf->prevInScaf = 0;
				cnt = getCntBetween ( prevCtg, thisCtg );
				cnt->nextInScaf = NULL;
			}

			prevCtg = * ( unsigned int * ) darrayGet ( tempArray, index2 - 1 );
			thisCtg = * ( unsigned int * ) darrayGet ( tempArray, index2 );
			cnt = getCntBetween ( prevCtg, thisCtg );

			if ( cnt->nextInScaf )
			{
				nextCtg = cnt->nextInScaf->contigID;
				cnt->nextInScaf->prevInScaf = 0;
				cnt = getCntBetween ( getTwinCtg ( nextCtg ), getTwinCtg ( thisCtg ) );
				cnt->nextInScaf = NULL;
			}

			prevCtg = * ( unsigned int * ) darrayGet ( tempArray, index1 );

			for ( t = index1 + 1; t <= index2; t++ )
			{
				thisCtg = * ( unsigned int * ) darrayGet ( tempArray, t );
				cnt = getCntBetween ( prevCtg, thisCtg );
				cnt->mask = 1;
				cnt->nextInScaf = NULL;
				cnt->prevInScaf = 0;
				cnt = getCntBetween ( getTwinCtg ( thisCtg ), getTwinCtg ( prevCtg ) );
				cnt->mask = 1;
				cnt->nextInScaf = NULL;
				cnt->prevInScaf = 0;
				prevCtg = thisCtg;
			}

			counter++;
		}
	}

	freeStack ( scafStack1 );
	freeStack ( scafStack2 );
	fprintf ( stderr, "Report from checkScaf: %d scaffold segments broken.\n", counter );
}

/*************************************************
Function:
    checkSimple
Description:
    Checks whether there is a contig appearing more than once
    in array.
Input:
    1. ctgArray:        array of contigs
    2. count:           contig number
Output:
    None.
Return:
    1 if no contig appeared more than once.
*************************************************/
static boolean checkSimple ( DARRAY * ctgArray, int count )
{
	int i;
	unsigned int ctg;

	for ( i = 0; i < count; i++ )
	{
		ctg = * ( unsigned int * ) darrayGet ( ctgArray, i );
		contig_array[ctg].flag = 0;
		contig_array[getTwinCtg ( ctg )].flag = 0;
	}

	for ( i = 0; i < count; i++ )
	{
		ctg = * ( unsigned int * ) darrayGet ( ctgArray, i );

		if ( contig_array[ctg].flag )
			{ return 0; }

		contig_array[ctg].flag = 1;
		contig_array[getTwinCtg ( ctg )].flag = 1;
	}

	return 1;
}

/*************************************************
Function:
    checkCircle
Description:
    Detects and masks contigs having circle connections.
    A --> B, B --> A.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void checkCircle()
{
	unsigned int i, ctg;
	CONNECT * cn_temp1;
	int counter = 0;

	for ( i = 1; i <= num_ctg; i++ )
	{
		cn_temp1 = contig_array[i].downwardConnect;

		while ( cn_temp1 )
		{
			if ( cn_temp1->weak || cn_temp1->deleted )
			{
				cn_temp1 = cn_temp1->next;
				continue;
			}

			ctg = cn_temp1->contigID;

			if ( checkConnect ( ctg, i ) )
			{
				counter++;
				maskContig ( i, 1 );
				maskContig ( ctg, 1 );
			}

			cn_temp1 = cn_temp1->next;
		}
	}

	fprintf ( stderr, "%d circles removed.\n", counter );
}
