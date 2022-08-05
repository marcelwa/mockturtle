/****************************************************************************
/*                                                                         */
/*                            QUICKCROSS                                   */
/*                                                                         */
/*                                                                         */
/*  Written by: Michael Haythorpe and Alex Newcombe                        */
/*  Date:       April 26, 2018                                             */
/*  Version:    1.0                                                        */
/*  Email:      michael.haythorpe@flinders.edu.au                          */
/*                                                                         */
/*                                                                         */
/*  Permission is granted for academic research use.  For other uses,      */
/*  contact the authors for licensing options.                             */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/*                                                                         */
/*                        LIST OF SUBROUTINES                              */
/*                                                                         */
/*  Runner Subroutines:                                                    */
/*                                                                         */
/*  main                     - execute program, read in graph data         */
/*  Print_Help               - print out command line argument options     */
/*  Create_Graph             - create a graph using the -c option          */
/*  Read_Graph               - read graph from file (edgelist format)      */
/*  Create_Graph6            - read graph6 code using the -g option        */
/*  Biconnected_Runner       - begin preprocessing and execute algorithm   */
/*  Find_One_Connected_Verts - find 1-connected vertices in preprocessing  */
/*                                                                         */
/*                                                                         */
/*  Initial Embedding Subroutines:                                         */
/*                                                                         */
/*  Embed_Circ               - use circle embedding                        */
/*  Embed_Planar_Subgraph    - find planar subgraph, then add vertices     */
/*                             using stripped down version of algorithm    */
/*  Planarise_Sixidx_PE      - stripped down version of Planarise_Sixidx   */
/*  Get_Distances_PE         - stripped down version of Get_Distances      */
/*  New_Clabel_PE            - stripped down version of New_Clabel         */
/*  Subd_Try_PE              - stripped down version of Subd_Try           */
/*  KKSL                     - Kamada Kawai spring layout algiruthm        */
/*  Get_Spring_Data          - convert vertex coordinates to an embedding  */
/*  UDFS                     - DFS search to find biconnected components   */
/*  UDFS5                    - DFS search used for embed scheme 5          */
/*  UDFS4                    - DFS search used for embed scheme 4          */
/*  DFSrelabel4              - relabel vertices for DFS                    */
/*                                                                         */
/*                                                                         */
/*  Planarity Check/Embedding Algorithms (for embed scheme 4):             */
/*                                                                         */
/*  Segment_Bipartite        - check if a segment is bipartite             */
/*  Strongly_Planar          - check if a subgraph is strongly planar      */
/*  Embedding_Given_Alpha    - construct an embedding from data from       */
/*                             Tarjan Hopcroft planarity check algorithm   */
/*  Embedding_Given_Crossings- construct an embedding from a crossing list */
/*                             and run Tarjan Hopcroft planarity check to  */
/*                             confirm that it's valid                     */
/*                                                                         */
/*                                                                         */
/*  Main Iteration Subroutines:                                            */
/*                                                                         */
/*  Quick_Cross_Main_Loop    - perform individual iterations of algorithm  */
/*  Planarise_Sixidx         - produce planarised graph                    */
/*  Faces                    - identify faces, and produce dual graph      */
/*  Prep_Vars                - set up variables for inner iteration        */
/*  Merge_Dual               - deleting vertices in original graph means   */
/*                             merging faces, hence merge dual vertices    */
/*  Get_Distances            - find relevant shortest paths in dual graph  */
/*  Sp_Bigface               - shortest path code for bigface scheme       */
/*  Find_Merge_Faces         - process shortest path through merged faces  */
/*  New_Clabel               - find new crossing labels after iteration    */
/*  Paths_Tree               - if multiple edges in new path cross same    */
/*                             edge, produce a tree to indicate orer       */
/*  Traverse_Paths_Tree      - interpret tree found in Paths_Tree          */
/*  New_Crossing_Orientation - determine orientation of a new crossing     */
/*  Subd_Try                 - see if any subdivisions are needed, also    */
/*                             consider incident edge movements            */
/*  Undo_Sub                 - see if any old subdivisions can be removed  */
/*                                                                         */
/*                                                                         */
/*  Utility Subroutines:                                                   */
/*                                                                         */
/*  dec_to_bin               - used to process graph6 format               */
/*  dec_to_dec               - used to process graph6 format               */
/*  randperm                 - produces a random permutation of integers   */
/*  ascend_cmpfunc_double    - comparison function for qsort               */
/*  descend_cmpfunc_double   - comparison function for qsort               */
/*  ascend_cmpfunc           - comparison function for qsort               */
/*  copy_array               - copy array of integers                      */
/*  copy_array_bool          - copy array of bools                         */
/*                                                                         */
/*                                                                         */
/***************************************************************************/

#pragma once

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void dec_to_bin( int n, bool* b )
{
  int i;
  int j;

  for ( i = 0; i < 6; i++ )
  {
    int num = 1;
    for ( j = 5; j > i; j-- )
      num = 2 * num;
    if ( n >= num )
    {
      b[i] = true;
      n = n - num;
    }
    else
      b[i] = false;
  }
}

int dec_bin_dec( int* n, int num )
{
  int i;

  int N = 0;
  if ( num == 3 )
  {
    bool b1[6];
    bool b2[6];
    bool b3[6];
    dec_to_bin( n[0], b1 );
    dec_to_bin( n[1], b2 );
    dec_to_bin( n[2], b3 );
    int multiplier = 1;
    for ( i = 0; i < 6; i++ )
    {
      if ( b3[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b2[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b1[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
  }
  if ( num == 6 )
  {
    bool b1[6];
    bool b2[6];
    bool b3[6];
    bool b4[6];
    bool b5[6];
    bool b6[6];
    dec_to_bin( n[0], b1 );
    dec_to_bin( n[1], b2 );
    dec_to_bin( n[2], b3 );
    dec_to_bin( n[4], b4 );
    dec_to_bin( n[5], b5 );
    dec_to_bin( n[6], b6 );
    int multiplier = 1;
    for ( i = 0; i < 6; i++ )
    {
      if ( b6[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b5[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b4[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b3[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b2[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
    for ( i = 0; i < 6; i++ )
    {
      if ( b1[5 - i] )
        N = N + multiplier;
      multiplier = 2 * multiplier;
    }
  }

  return N;
}

void randperm( int** rp, int n )
{
  int i;
  for ( i = 0; i < n; i++ )
    ( *rp )[i] = ( i + 1 );

  for ( i = 0; i < n; i++ )
  {
    int j, t;
    j = rand() % ( n - i ) + i;
    t = ( *rp )[j];
    ( *rp )[j] = ( *rp )[i];
    ( *rp )[i] = t; // Swap i and j
  }
}

int ascend_cmpfunc_double( const void* a, const void* b )
{
  if ( *(double*)a == *(double*)b )
    return 0;
  if ( *(double*)a - *(double*)b > 0 )
    return 1;
  else
    return -1;
}

int descend_cmpfunc_double( const void* a, const void* b )
{
  if ( *(double*)a == *(double*)b )
    return 0;
  if ( *(double*)a - *(double*)b < 0 )
    return 1;
  else
    return -1;
}

int ascend_cmpfunc( const void* a, const void* b )
{
  if ( *(int*)a == *(int*)b )
    return 0;
  if ( *(int*)a - *(int*)b > 0 )
    return 1;
  else
    return -1;
}

void copy_array( int* a1, int** a2, int n )
{
  int i;
  *a2 = (int*)malloc( n * sizeof( int ) );
  for ( i = 0; i < n; i++ )
    ( *a2 )[i] = a1[i];
}

void copy_array_bool( bool* a1, bool** a2, int n )
{
  int i;
  *a2 = (bool*)malloc( n * sizeof( bool ) );
  for ( i = 0; i < n; i++ )
    ( *a2 )[i] = a1[i];
}

int Undo_Sub( int N, int M, int** Sixidx, int* Clabel, int* Dlabel, int** Cidx, int** incidente, int** deg, int N_before_Subd, int orig_N, int verbose )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Undo_Sub                                                               */
  /*                                                                         */
  /*  This subroutine considers any previous subdivided vertices to see if   */
  /*  they can now safely be removed. Subdivisions are only needed if their  */
  /*  removal would mean one edge crosses another edge multiple times.       */
  /*                                                                         */
  /*  Vertex labels are such that vertices resulting from subdivision all    */
  /*  have larger labels than original vertices. Only these need checking.   */
  /*  There is no need to check subdivided vertices made in the last step.   */
  /*                                                                         */
  /*  If an unnecessary subdivided vertex is found, it is removed and the    */
  /*  two incident edges are merged into one.                                */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;
  int total_removed = 0;

  bool* checked;
  if ( N_before_Subd > orig_N )
  {
    checked = (bool*)malloc( M * sizeof( bool ) );
    for ( i = 0; i < M; i++ )
      checked[i] = false;
  }
  else
  {
    // fprintf(stdout,"\n exiting early \n");
    return total_removed;
  }

  for ( i = N_before_Subd; i >= orig_N + 1; i-- )
  {
    // First check to see if subdivision is still required. Vertex will be degree 2, get
    // both edges and see if they are both crossed by a common edge.

    bool asmall = true;
    bool bsmall = true;

    int e1 = ( *incidente )[i - 1];
    int a = ( *Sixidx )[e1 - 1];
    int e2 = ( *Sixidx )[e1 - 1 + 4 * M];
    if ( a == i )
    {
      asmall = false;
      a = ( *Sixidx )[e1 - 1 + M];
      e2 = ( *Sixidx )[e1 - 1 + 2 * M];
    }

    int b = ( *Sixidx )[e2 - 1];
    if ( b == i )
    {
      bsmall = false;
      b = ( *Sixidx )[e2 - 1 + M];
    }

    for ( j = ( *Cidx )[e1 - 1]; j < ( *Cidx )[e1]; j++ )
      checked[Clabel[j] - 1] = true;

    bool remove = true;

    for ( j = ( *Cidx )[e2 - 1]; j < ( *Cidx )[e2]; j++ )
    {
      if ( checked[Clabel[j] - 1] )
      {
        remove = false;
        break;
      }
      if ( Clabel[j] == e1 )
      {
        remove = false;
        break;
      }
    }
    int orige1 = e1;

    if ( remove )
    {
      // Subdivision is no longer required. Go through process of removing
      // the vertex and merging the two edges into one.

      total_removed++;

      if ( a > b )
      {
        int temp = a;
        a = b;
        b = temp;
        bool tempsmall;
        tempsmall = asmall;
        asmall = bsmall;
        bsmall = tempsmall;
        temp = e1;
        e1 = e2;
        e2 = temp;
      }

      // Now we have a < b, and asmall = 1 iff a < v, bsmall = 1 iff b < v

      // Need to:
      // Check if edge is av or va.
      //
      // If av and other edge is bv:
      //        Change edge to ab, need to update Sixidx.
      //        Delete edge bv, moving Clabel and Dlabel into edge ab. Orders reversed.
      //        Sixidx(ab,5:6) becomes [Sixidx(bv,3) Sixidx(bv,4)]
      //
      // If av and other edge is vb:
      //        Change edge to ab, need to update Sixidx.
      //        Delete edge vb, moving Clabel and Dlabel into edge ab. Order stays the same.
      //        Sixidx(ab,5:6) becomes [Sixidx(vb,5) Sixidx(vb,6)]
      //
      // If va: Change edge to ab, need to update Sixidx, and turn Clabel around, and opposite Dlabels
      //        Delete edge vb, moving Clabel and Dlabel into edge ab. Order stays the same.
      //        Sixidx(ab,3:6) becomes [Sixidx(ab,5) Sixidx(ab,6) Sixidx(vb,5) Sixidx(vb,6)]
      // Need to modify Cidx from a to b, and resize
      // Incidente for b turns into edge ab
      // Remove v (means resizing incidente, deg)

      int e1c = ( *Cidx )[e1] - ( *Cidx )[e1 - 1];
      int e2c = ( *Cidx )[e2] - ( *Cidx )[e2 - 1];
      if ( asmall )
      {
        if ( bsmall )
        {
          ( *Sixidx )[e1 - 1 + M] = b;
          ( *Sixidx )[e1 - 1 + 4 * M] = ( *Sixidx )[e2 - 1 + 2 * M];
          ( *Sixidx )[e1 - 1 + 5 * M] = ( *Sixidx )[e2 - 1 + 3 * M];
        }
        else
        {
          ( *Sixidx )[e1 - 1 + M] = b;
          ( *Sixidx )[e1 - 1 + 4 * M] = ( *Sixidx )[e2 - 1 + 4 * M];
          ( *Sixidx )[e1 - 1 + 5 * M] = ( *Sixidx )[e2 - 1 + 5 * M];
        }
      }
      else
      {
        ( *Sixidx )[e1 - 1] = a;
        ( *Sixidx )[e1 - 1 + M] = b;
        ( *Sixidx )[e1 - 1 + 2 * M] = ( *Sixidx )[e1 - 1 + 4 * M];
        ( *Sixidx )[e1 - 1 + 3 * M] = ( *Sixidx )[e1 - 1 + 5 * M];
        ( *Sixidx )[e1 - 1 + 4 * M] = ( *Sixidx )[e2 - 1 + 4 * M];
        ( *Sixidx )[e1 - 1 + 5 * M] = ( *Sixidx )[e2 - 1 + 5 * M];

        int* temp_Clabel_e1 = (int*)malloc( e1c * sizeof( int ) );
        int* temp_Dlabel_e1 = (int*)malloc( e1c * sizeof( int ) );
        for ( j = ( *Cidx )[e1 - 1]; j < ( *Cidx )[e1]; j++ )
        {
          temp_Clabel_e1[j - ( *Cidx )[e1 - 1]] = Clabel[j];
          temp_Dlabel_e1[j - ( *Cidx )[e1 - 1]] = Dlabel[j];
        }

        for ( j = 0; j < e1c; j++ )
        {
          Clabel[( *Cidx )[e1 - 1] + j] = temp_Clabel_e1[e1c - j - 1];
          Dlabel[( *Cidx )[e1 - 1] + j] = -temp_Dlabel_e1[e1c - j - 1];
          int oe = Clabel[( *Cidx )[e1 - 1] + j];
          for ( k = ( *Cidx )[oe - 1]; k < ( *Cidx )[oe]; k++ )
          {
            if ( Clabel[k] == e1 )
            {
              Dlabel[k] = -Dlabel[k];
              break;
            }
          }
        }

        free( temp_Clabel_e1 );
        free( temp_Dlabel_e1 );
      }

      int* temp_Clabel_e2 = (int*)malloc( e2c * sizeof( int ) );
      int* temp_Dlabel_e2 = (int*)malloc( e2c * sizeof( int ) );
      for ( j = ( *Cidx )[e2 - 1]; j < ( *Cidx )[e2]; j++ )
      {
        temp_Clabel_e2[j - ( *Cidx )[e2 - 1]] = Clabel[j];
        temp_Dlabel_e2[j - ( *Cidx )[e2 - 1]] = Dlabel[j];
      }
      if ( e1 < e2 )
      {
        // Everything from e1+1 to e2-1 moves forward e2c places
        for ( j = ( *Cidx )[e2 - 1] - 1; j >= ( *Cidx )[e1]; j-- )
        {
          Clabel[j + e2c] = Clabel[j];
          Dlabel[j + e2c] = Dlabel[j];
        }

        // Every Cidx from e1 to e2-1 increases by e2c, every Cidx from e2 to M is set to following one
        for ( j = e1; j <= e2 - 1; j++ )
          ( *Cidx )[j] = ( *Cidx )[j] + e2c;
        for ( j = e2; j < M; j++ )
          ( *Cidx )[j] = ( *Cidx )[j + 1];
        *Cidx = (int*)realloc( *Cidx, M * sizeof( int ) );
      }
      else
      {
        // Everything from e2+1 to e1 moves backwards e2c places
        for ( j = ( *Cidx )[e2]; j < ( *Cidx )[e1]; j++ )
        {
          Clabel[j - e2c] = Clabel[j];
          Dlabel[j - e2c] = Dlabel[j];
        }

        // Every Cidx from e2 to e1-1 is equal to next one minus e2c, every Cidx from e1 to M is set to following one
        for ( j = e2; j <= e1 - 2; j++ )
          ( *Cidx )[j] = ( *Cidx )[j + 1] - e2c;
        for ( j = e1 - 1; j < M; j++ )
          ( *Cidx )[j] = ( *Cidx )[j + 1];

        *Cidx = (int*)realloc( *Cidx, M * sizeof( int ) );
      }
      int offset = 0;
      if ( e1 > e2 )
        offset = 1;

      if ( bsmall )
      {

        for ( j = 0; j < e2c; j++ )
        {
          Clabel[( *Cidx )[e1 - 1 - offset] + e1c + j] = temp_Clabel_e2[e2c - j - 1];
          Dlabel[( *Cidx )[e1 - 1 - offset] + e1c + j] = -temp_Dlabel_e2[e2c - j - 1];
        }
      }
      else
      {
        for ( j = 0; j < e2c; j++ )
        {
          Clabel[( *Cidx )[e1 - 1 - offset] + e1c + j] = temp_Clabel_e2[j];
          Dlabel[( *Cidx )[e1 - 1 - offset] + e1c + j] = temp_Dlabel_e2[j];
        }
      }

      free( temp_Clabel_e2 );
      free( temp_Dlabel_e2 );
      ( *incidente )[b - 1] = e1;

      // Now resizing incidente and deg
      for ( j = i; j <= N - 1; j++ )
      {
        ( *incidente )[j - 1] = ( *incidente )[j];
        ( *deg )[j - 1] = ( *deg )[j];
      }

      *incidente = (int*)realloc( *incidente, ( N - 1 ) * sizeof( int ) );
      *deg = (int*)realloc( *deg, ( N - 1 ) * sizeof( int ) );

      // Remove row e2 from Sixidx
      int increment = 0;
      for ( j = e2 - 1; j < 6 * M - 6; j++ )
      {
        if ( j % M == ( e2 - 1 - increment ) % M )
          increment++;
        ( *Sixidx )[j] = ( *Sixidx )[j + increment];
      }

      *Sixidx = (int*)realloc( *Sixidx, 6 * ( M - 1 ) * sizeof( int ) );

      N = N - 1;
      M = M - 1;

      // Any vertex labels in Sixidx larger than i need to be subtracted by 1
      for ( j = 0; j < 2 * M; j++ )
        if ( ( *Sixidx )[j] > i )
          ( *Sixidx )[j] = ( *Sixidx )[j] - 1;

      // Finally, anything in Sixidx, Clabel or incidente equal to e2 changes to e1. Anything bigger than e2 gets subtracted by 1.
      for ( j = 2 * M; j < 6 * M; j++ )
      {
        if ( ( *Sixidx )[j] == e2 )
          ( *Sixidx )[j] = e1;
        if ( ( *Sixidx )[j] > e2 )
          ( *Sixidx )[j] = ( *Sixidx )[j] - 1;
      }

      for ( j = 0; j < ( *Cidx )[M]; j++ )
      {
        if ( Clabel[j] == e2 )
        {
          Clabel[j] = e1;
          if ( bsmall )
            Dlabel[j] = -Dlabel[j];
        }
        if ( Clabel[j] > e2 )
          Clabel[j] = Clabel[j] - 1;
      }

      for ( j = 0; j < N; j++ )
        if ( ( *incidente )[j] > e2 )
          ( *incidente )[j] = ( *incidente )[j] - 1;

      if ( verbose > 1 )
        fprintf( stdout, "Removing subdivided vertex %d, now N = %d\n", i, N );
    }

    if ( i > orig_N + 1 )
      for ( j = 0; j < M; j++ )
        checked[j] = false;
  }

  free( checked );

  return total_removed;
}

void Subd_Try( int current_crossings, int** Sixidx, int* Clabel, int* Dlabel, int** Cidx, int* pathedge, int N, int M, int* origedgelabel, int** incidente, int** deg, int* pdist, int* fi, int* fin, int vm, int d, int maxpdist, int* pidx, int* outputs, int orig_N, bool only_incident_check, int verbose )
{ // IN ADVANCE:

  /***************************************************************************/
  /*                                                                         */
  /*  Subd_Try                                                               */
  /*                                                                         */
  /*  This subroutine serves two purposes. First it checks if any edge (a,b) */
  /*  immediately crosses an edge (a,d). This will not happen commonly, but  */
  /*  can happen after another vertex is moved out of the way. We call this  */
  /*  an "incident edge crossing". If it happens, edge (a,b) is relocated so */
  /*  as to not cross (a,d) and the number of crossings is updated.          */
  /*                                                                         */
  /*  The second purpose is to check if any subdivisions are required. They  */
  /*  are required if any edge crosses another edge multiple times. This     */
  /*  causes issues later, and so we avoid it by subdividing the edge each   */
  /*  it crosses the same edge.                                              */
  /*                                                                         */
  /*  Afterwards, Undo_Sub is called to see if any old subdivisions can now  */
  /*  be removed. If they are removed, Subd_Try is run again since it is     */
  /*  possible that new incident edge crossings have been created. The check */
  /*  for subdivisions is skipped in this situation.                         */
  /*                                                                         */
  /***************************************************************************/

  int N_before_Subd = N;

  int curre;
  int othere;
  int i;
  int j;
  int k;
  int tt;

  // First, check for incident edge crossings.

  bool consider_incidente = true;
  while ( consider_incidente )
  {
    consider_incidente = false;
    for ( i = 0; i < N; i++ )
    {
      curre = ( *incidente )[i];
      for ( j = 0; j < ( *deg )[i]; j++ )
      {
        if ( ( *Cidx )[curre] - ( *Cidx )[curre - 1] > 0 )
        {
          if ( ( *Sixidx )[curre - 1] == i + 1 )
          {
            if ( Clabel[( *Cidx )[curre - 1]] == ( *Sixidx )[( 2 * M ) - 1 + curre] )
            {
              othere = ( *Sixidx )[( 2 * M ) - 1 + curre];
              if ( ( *Sixidx )[othere - 1] == i + 1 && Clabel[( *Cidx )[othere - 1]] == curre )
              {
                consider_incidente = true;
                for ( k = ( *Cidx )[curre - 1]; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = curre + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }
                for ( k = ( *Cidx )[othere - 1]; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = othere + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }

                if ( ( *deg )[i] > 2 )
                {
                  if ( ( *Sixidx )[( *Sixidx )[( 2 * M ) - 1 + othere] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                  }
                  else
                  {
                    ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                  }
                  if ( ( *Sixidx )[( *Sixidx )[( 3 * M ) - 1 + curre] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                  }
                  else
                  {
                    ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                  }
                  tt = ( *Sixidx )[( 3 * M ) - 1 + curre];
                  ( *Sixidx )[( 2 * M ) - 1 + curre] = ( *Sixidx )[( 2 * M ) - 1 + othere];
                  ( *Sixidx )[( 3 * M ) - 1 + curre] = othere;
                  ( *Sixidx )[( 2 * M ) - 1 + othere] = curre;
                  ( *Sixidx )[( 3 * M ) - 1 + othere] = tt;
                }

                current_crossings = current_crossings - 1;
              }
              else if ( ( *Sixidx )[M - 1 + othere] == i + 1 && Clabel[( *Cidx )[othere] - 1] == curre )
              {
                consider_incidente = true;
                for ( k = ( *Cidx )[curre - 1]; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = curre + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }
                for ( k = ( *Cidx )[othere] - 1; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = othere + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }

                if ( ( *deg )[i] > 2 )
                {
                  if ( ( *Sixidx )[( *Sixidx )[( 4 * M ) - 1 + othere] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 3 * M ) + ( *Sixidx )[( 4 * M ) - 1 + othere] - 1] = curre;
                  }
                  else
                  {
                    ( *Sixidx )[( 5 * M ) + ( *Sixidx )[( 4 * M ) - 1 + othere] - 1] = curre;
                  }
                  if ( ( *Sixidx )[( *Sixidx )[( 3 * M ) - 1 + curre] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                  }
                  else
                  {
                    ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                  }
                  tt = ( *Sixidx )[( 3 * M ) - 1 + curre];
                  ( *Sixidx )[( 2 * M ) - 1 + curre] = ( *Sixidx )[( 4 * M ) - 1 + othere];
                  ( *Sixidx )[( 3 * M ) - 1 + curre] = othere;
                  ( *Sixidx )[( 4 * M ) - 1 + othere] = curre;
                  ( *Sixidx )[( 5 * M ) - 1 + othere] = tt;
                }
                current_crossings = current_crossings - 1;
              }
            }
          }
          else
          {

            if ( Clabel[( *Cidx )[curre] - 1] == ( *Sixidx )[( 4 * M ) - 1 + curre] )
            {
              othere = ( *Sixidx )[( 4 * M ) - 1 + curre];
              if ( ( *Sixidx )[othere - 1] == i + 1 && Clabel[( *Cidx )[othere - 1]] == curre )
              {
                consider_incidente = true;
                for ( k = ( *Cidx )[curre] - 1; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = curre + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }
                for ( k = ( *Cidx )[othere - 1]; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = othere + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }

                if ( ( *deg )[i] > 2 )
                {
                  if ( ( *Sixidx )[( *Sixidx )[( 2 * M ) - 1 + othere] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                  }
                  else
                  {
                    ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                  }
                  if ( ( *Sixidx )[( *Sixidx )[( 5 * M ) - 1 + curre] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                  }
                  else
                  {
                    ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                  }
                  tt = ( *Sixidx )[( 5 * M ) - 1 + curre];
                  ( *Sixidx )[( 4 * M ) - 1 + curre] = ( *Sixidx )[( 2 * M ) - 1 + othere];
                  ( *Sixidx )[( 5 * M ) - 1 + curre] = othere;
                  ( *Sixidx )[( 2 * M ) - 1 + othere] = curre;
                  ( *Sixidx )[( 3 * M ) - 1 + othere] = tt;
                }
                current_crossings = current_crossings - 1;
              }
              else if ( ( *Sixidx )[M - 1 + othere] == i + 1 && Clabel[( *Cidx )[othere] - 1] == curre )
              {
                consider_incidente = true;
                for ( k = ( *Cidx )[curre] - 1; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = curre + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }
                for ( k = ( *Cidx )[othere] - 1; k < ( *Cidx )[M] - 1; k++ )
                {
                  Clabel[k] = Clabel[k + 1];
                  Dlabel[k] = Dlabel[k + 1];
                }
                for ( k = othere + 1; k < M + 2; k++ )
                {
                  ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                }

                if ( ( *deg )[i] > 2 )
                {
                  if ( ( *Sixidx )[( *Sixidx )[( 4 * M ) - 1 + othere] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + othere]] = curre;
                  }
                  else
                  {
                    ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + othere]] = curre;
                  }
                  if ( ( *Sixidx )[( *Sixidx )[( 5 * M ) - 1 + curre] - 1] == i + 1 )
                  {
                    ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                  }
                  else
                  {
                    ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                  }
                  tt = ( *Sixidx )[( 5 * M ) - 1 + curre];
                  ( *Sixidx )[( 4 * M ) - 1 + curre] = ( *Sixidx )[( 4 * M ) - 1 + othere];
                  ( *Sixidx )[( 5 * M ) - 1 + curre] = othere;
                  ( *Sixidx )[( 4 * M ) - 1 + othere] = curre;
                  ( *Sixidx )[( 5 * M ) - 1 + othere] = tt;
                }
                current_crossings = current_crossings - 1;
              }
            }
          }
        }
        if ( ( *Sixidx )[curre - 1] == i + 1 )
        {
          curre = ( *Sixidx )[( 2 * M ) + curre - 1];
        }
        else
        {
          curre = ( *Sixidx )[( 4 * M ) - 1 + curre];
        }
      }
    }
  }

  if ( only_incident_check == false )
  {

    // Now check for subdivisions. These will occur if an edge crosses the same edge multiple times.

    int restartloop = 1;
    int* pcnt = (int*)malloc( d * sizeof( int ) );
    int oldN = N;
    int oldM = M;
    int ci;
    int qi;
    int qi1;
    int cj;
    int turnchk1;
    int chksum = 0;
    int r;
    int* cc = (int*)malloc( maxpdist * sizeof( int ) );
    int v2;
    int* newc1 = (int*)malloc( maxpdist * sizeof( int ) );
    int* newc2 = (int*)malloc( maxpdist * sizeof( int ) );
    int* newd1 = (int*)malloc( maxpdist * sizeof( int ) );
    int* newd2 = (int*)malloc( maxpdist * sizeof( int ) );
    int ncnt1;
    int ncnt2;
    int des;
    int cnt3;
    int cnt4;
    int qc;
    int qf;
    int* pe = (int*)malloc( maxpdist * sizeof( int ) );
    int* ze = (int*)malloc( maxpdist * sizeof( int ) );
    int* de = (int*)malloc( maxpdist * sizeof( int ) );
    int erg;
    int qrg;
    int cntr;
    int cntr2;
    int* qe = (int*)malloc( maxpdist * sizeof( int ) );
    int q;
    int* Mid = (int*)malloc( d * maxpdist * sizeof( int ) ); // size of vector here was causing issues? - originally maxpdist, but not long enough?? i think sumpdist is the correct size
    int* oe = (int*)malloc( maxpdist * sizeof( int ) );

    while ( restartloop == 1 )
    {
      restartloop = 0;

      for ( i = 1; i < d + M - oldM + 1; i++ )
      {
        if ( i <= d )
        {
          ci = fi[i - 1];
        }
        else
        {
          ci = oldM + i - d;
        }
        qi = ( *Cidx )[ci - 1];
        qi1 = ( *Cidx )[ci];
        if ( ( *Sixidx )[M + ci - 1] > oldN )
        {
          turnchk1 = 1;
        }
        else
        {
          turnchk1 = 0;
        }
        for ( j = 0; j < qi1 - qi; j++ )
        {
          if ( turnchk1 == 1 )
          {
            cj = qi1 - j - 1;
          }
          else
          {
            cj = qi + j;
          }
          chksum = 0;
          for ( k = qi; k < qi1; k++ )
          {
            if ( Clabel[k] == Clabel[cj] )
            {
              chksum = chksum + 1;
            }
          }
          if ( chksum > 1 )
          {
            // then subdvision incoming so reallocate (*Sixidx), (*Cidx), (*incidente), (*deg)
            *Sixidx = (int*)realloc( *Sixidx, ( 6 * M + 6 ) * sizeof( int ) );
            for ( r = 0; r < 5; r++ )
            {
              for ( k = r * M; k < ( r + 1 ) * M; k++ )
              {
                ( *Sixidx )[6 * M + 6 - k - 2 - r] = ( *Sixidx )[6 * M - k - 1];
              }
            }
            for ( k = 0; k < 5; k++ )
            {
              ( *Sixidx )[( k + 1 ) * ( M + 1 ) - 1] = 0;
            }
            *incidente = (int*)realloc( *incidente, ( N + 1 ) * sizeof( int ) );
            *deg = (int*)realloc( *deg, ( N + 1 ) * sizeof( int ) );
            *Cidx = (int*)realloc( *Cidx, ( M + 2 ) * sizeof( int ) );
            for ( k = 0; k < maxpdist; k++ )
            { // reset cc for use again
              cc[k] = 0;
            }
            v2 = Clabel[cj];
            if ( verbose > 1 )
              fprintf( stdout, "Subdividing an edge: %d with %d\n", ci, v2 );

            N = N + 1;
            M = M + 1;
            cntr = 0;
            cntr2 = 0;
            ncnt1 = -1;
            ncnt2 = -1;

            if ( ( *Sixidx )[M - 1 + ci] <= oldN )
            {
              for ( k = qi; k < cj + 1; k++ )
              {
                ncnt1 = ncnt1 + 1;
                newc1[ncnt1] = Clabel[k];
                newd1[ncnt1] = Dlabel[k];
              }
              for ( k = 0; k < ( qi1 - cj - 1 ); k++ )
              {

                newc2[k] = Clabel[qi1 - 1 - k];
                newd2[k] = -Dlabel[qi1 - 1 - k];
              }
              ncnt2 = ( qi1 - cj - 2 );
            }
            else
            {
              for ( k = 0; k < ( qi1 - cj ); k++ )
              {
                ncnt1 = ncnt1 + 1;
                newc1[k] = Clabel[qi1 - 1 - k];
                newd1[k] = -Dlabel[qi1 - 1 - k];
              }
              for ( k = qi; k < cj; k++ )
              {

                ncnt2 = ncnt2 + 1;
                newc2[ncnt2] = Clabel[k];
                newd2[ncnt2] = Dlabel[k];
              }
            }
            cntr = -1;
            qc = ( *Cidx )[v2 - 1];
            for ( k = qc; k < ( *Cidx )[v2]; k++ )
            {
              if ( Clabel[k] == ci )
              {
                cntr = cntr + 1;
                cc[cntr] = k - qc;
              }
            }

            for ( k = 0; k < maxpdist; k++ )
            {
              pe[k] = 0;
              qe[k] = 0;
              ze[k] = 0;
            }

            cntr = -1;
            if ( ( *Sixidx )[M - 1 + ci] > oldN )
            {
              q = Mid[ci - oldM - 1];
              if ( ( *Sixidx )[ci - 1] == fin[q - 1] )
              {
                for ( k = pcnt[q - 1] + 1; k < pdist[q - 1] + 1; k++ )
                {
                  cntr = cntr + 1;
                  pe[cntr] = pathedge[pidx[q - 1] + k - 1];
                }
              }
              else
              {
                for ( k = 1; k < pcnt[q - 1]; k++ )
                {
                  cntr = cntr + 1;
                  pe[cntr] = pathedge[pidx[q - 1] + pcnt[q - 1] - k - 1];
                }
              }
            }
            else
            {
              for ( k = 0; k < d; k++ )
              {
                if ( fi[k] == ci )
                {
                  qf = k;
                  break;
                }
              }
              if ( ( *Sixidx )[ci - 1] == vm )
              {
                for ( k = pidx[qf]; k < pidx[qf + 1]; k++ )
                {
                  cntr = cntr + 1;
                  pe[cntr] = pathedge[k];
                }
              }
              else
              {
                for ( k = 0; k < pdist[qf]; k++ )
                {
                  cntr = cntr + 1;
                  pe[cntr] = pathedge[pidx[qf + 1] - 1 - k];
                }
              }
            }

            cntr2 = -1;
            for ( k = 0; k < cntr + 1; k++ )
            {
              oe[k] = origedgelabel[pe[k] - 1];
              if ( oe[k] == v2 )
              {

                cntr2 = cntr2 + 1;
                ze[cntr2] = k + 1;
                de[cntr2] = pe[k];
              }
            }

            des = de[0];

            qsort( de, cntr2 + 1, sizeof( int ), ascend_cmpfunc );

            if ( ( *Sixidx )[M - 1 + ci] > oldN )
            {
              Mid[M - oldM - 1] = q;
              if ( ( *Sixidx )[ci - 1] == fin[q - 1] )
              {
                pcnt[q - 1] = pcnt[q - 1] + ze[0];
              }
              else
              {
                pcnt[q - 1] = pcnt[q - 1] - ze[0];
              }
            }
            else
            {
              Mid[M - oldM - 1] = qf + 1;
              if ( ( *Sixidx )[ci - 1] == vm )
              {
                pcnt[qf] = ze[0];
              }
              else
              {
                pcnt[qf] = pdist[qf] + 1 - ze[0];
              }
            }

            bool lessN = false;
            if ( ( *Sixidx )[M - 1 + ci] <= oldN )
              lessN = true;

            for ( k = 0; k < cntr2 + 1; k++ )
            {
              if ( de[k] != des )
              {
                Clabel[qc + cc[k]] = M;
                if ( lessN )
                {
                  Dlabel[qc + cc[k]] = -Dlabel[qc + cc[k]];
                }
              }
              else
              {
                if ( lessN == false )
                {
                  Dlabel[qc + cc[k]] = -Dlabel[qc + cc[k]];
                }
              }
            }

            if ( ( *Sixidx )[M - 1 + ci] <= oldN )
            {
              ( *Sixidx )[M - 1] = ( *Sixidx )[M - 1 + ci];
              ( *Sixidx )[( 2 * M ) - 1] = N;
              ( *Sixidx )[( 3 * M ) - 1] = ( *Sixidx )[( 4 * M ) - 1 + ci];
              ( *Sixidx )[( 4 * M ) - 1] = ( *Sixidx )[( 5 * M ) - 1 + ci];
              ( *Sixidx )[( 5 * M ) - 1] = ci;
              ( *Sixidx )[( 6 * M ) - 1] = ci;
              if ( ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + ci]] == ci )
              {
                ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + ci]] = M;
              }
              else
              {
                ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + ci]] = M;
              }
              if ( ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + ci]] == ci )
              {
                ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + ci]] = M;
              }
              else
              {
                ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + ci]] = M;
              }
              ( *Sixidx )[M - 1 + ci] = N;
              ( *Sixidx )[( 4 * M ) - 1 + ci] = M;
              ( *Sixidx )[( 5 * M ) - 1 + ci] = M;
            }
            else
            {
              ( *Sixidx )[M - 1] = ( *Sixidx )[ci - 1];
              ( *Sixidx )[( 2 * M ) - 1] = N;
              ( *Sixidx )[( 3 * M ) - 1] = ( *Sixidx )[( 2 * M ) - 1 + ci];
              ( *Sixidx )[( 4 * M ) - 1] = ( *Sixidx )[( 3 * M ) - 1 + ci];
              ( *Sixidx )[( 5 * M ) - 1] = ci;
              ( *Sixidx )[( 6 * M ) - 1] = ci;
              if ( ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + ci]] == ci )
              {
                ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + ci]] = M;
              }
              else
              {
                ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + ci]] = M;
              }
              if ( ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + ci]] == ci )
              {
                ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + ci]] = M;
              }
              else
              {
                ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + ci]] = M;
              }
              ( *Sixidx )[ci - 1] = ( *Sixidx )[M - 1 + ci];
              ( *Sixidx )[M - 1 + ci] = N;
              ( *Sixidx )[( 2 * M ) - 1 + ci] = ( *Sixidx )[( 4 * M ) - 1 + ci];
              ( *Sixidx )[( 3 * M ) - 1 + ci] = ( *Sixidx )[( 4 * M ) - 1 + ci];
              ( *Sixidx )[( 4 * M ) - 1 + ci] = M;
              ( *Sixidx )[( 5 * M ) - 1 + ci] = M;
            }

            for ( k = 0; k < ncnt2 + 1; k++ )
            {
              if ( newc2[k] != v2 )
              {
                erg = newc2[k];
                for ( r = ( *Cidx )[erg - 1]; r < ( *Cidx )[erg]; r++ )
                {
                  if ( Clabel[r] == ci )
                  {

                    Clabel[r] = M;
                    if ( lessN )
                      Dlabel[r] = -Dlabel[r];
                  }
                }
              }
            }
            if ( lessN == false )
            {
              for ( k = 0; k < ncnt1 + 1; k++ )
              {
                if ( newc1[k] != v2 )
                {
                  erg = newc1[k];
                  for ( r = ( *Cidx )[erg - 1]; r < ( *Cidx )[erg]; r++ )
                  {
                    if ( Clabel[r] == ci )
                    {

                      Dlabel[r] = -Dlabel[r];
                    }
                  }
                }
              }
            }
            cnt3 = -1;
            cnt4 = -1;
            for ( k = 0; k < ( *Cidx )[M - 1]; k++ )
            {
              if ( k >= qi && k <= qi + ncnt1 )
              { // these need checked
                cnt3 = cnt3 + 1;
                Clabel[k] = newc1[cnt3];
                Dlabel[k] = newd1[cnt3];
              }

              if ( k > qi + ncnt1 && k < ( *Cidx )[M - 1] - ncnt2 - 1 )
              {
                Clabel[k] = Clabel[k + ncnt2 + 1];
                Dlabel[k] = Dlabel[k + ncnt2 + 1];
              }
              if ( k >= ( *Cidx )[M - 1] - ncnt2 - 1 )
              {
                cnt4 = cnt4 + 1;
                Clabel[k] = newc2[cnt4];
                Dlabel[k] = newd2[cnt4];
              }
            }
            for ( k = ci; k < M; k++ )
            {
              ( *Cidx )[k] = ( *Cidx )[k] - cnt4 - 1;
            }
            ( *Cidx )[M] = ( *Cidx )[M - 1] + cnt4 + 1;
            ( *incidente )[N - 1] = M;
            ( *incidente )[( *Sixidx )[M - 1] - 1] = M;
            ( *deg )[N - 1] = 2;

            restartloop = 1;
            break;
          }
        }
        if ( restartloop == 1 )
        {
          break;
        }
      }
    }
    // fprintf(stdout,"cend000 \n");
    free( pcnt );
    free( cc );
    free( newc1 );
    free( newc2 );
    free( newd1 );
    free( newd2 );
    free( pe );
    free( ze );
    free( de );
    free( qe );
    free( Mid );
    free( oe );
    // fprintf(stdout,"cendxxx\n");

    // Next, see if any old subdivisions can be removed.
    // fprintf(stdout,"cend \n");

    int total_removed = Undo_Sub( N, M, Sixidx, Clabel, Dlabel, Cidx, incidente, deg, N_before_Subd, orig_N, verbose );

    if ( total_removed > 0 )
    {
      // If any old subdivisions are removed, check for incident edge crossings again.

      N = N - total_removed;
      M = M - total_removed;
      int newoutputs[3];
      Subd_Try( current_crossings, Sixidx, Clabel, Dlabel, Cidx, pathedge, N, M, origedgelabel, incidente, deg, pdist, fi, fin, vm, d, maxpdist, pidx, newoutputs, orig_N, true, verbose );
      N = newoutputs[0];
      M = newoutputs[1];
      current_crossings = newoutputs[2];
    }
  }

  outputs[0] = N;
  outputs[1] = M;
  outputs[2] = current_crossings;
}

int New_Crossing_Orientation( int Md, int* Sixidx_plan, int facesize, int* F, int* Fidx, int* origedgeorder, int ne, int vm, int vmN, int face, int org, int* elist, int* eidx )
{

  /***************************************************************************/
  /*                                                                         */
  /*  New_Crossing_Orientation                                               */
  /*                                                                         */
  /*  This subroutine determines the orientation of a new crossing.          */
  /*                                                                         */
  /***************************************************************************/

  int q = eidx[org] - eidx[org - 1];
  int f = -1;
  int i;
  int ori;
  for ( i = 0; i < facesize; i++ )
  {
    if ( F[Fidx[face - 1] + i] == ne )
    {
      f = i;
      break;
    }
  }
  int mf;
  if ( f > 0 )
  {
    mf = f - 1;
  }
  else
  {
    mf = facesize - 1;
  }
  int ck;
  int f1;
  if ( q > 1 )
  {
    f1 = origedgeorder[ne - 1];
    ck = 0;
    if ( f1 == 1 )
    {
      ck = 1;
    }
    else
    {
      if ( Sixidx_plan[elist[eidx[org - 1] + f1 - 2] - 1] == Sixidx_plan[ne - 1] || Sixidx_plan[Md + elist[eidx[org - 1] + f1 - 2] - 1] == Sixidx_plan[ne - 1] )
      {
        ck = 1;
      }
    }
    if ( ck == 1 )
    {
      if ( Sixidx_plan[F[Fidx[face - 1] + mf] - 1] == Sixidx_plan[ne - 1] || Sixidx_plan[Md + F[Fidx[face - 1] + mf] - 1] == Sixidx_plan[ne - 1] )
      {
        if ( vm < vmN )
        {
          ori = -1;
        }
        else
        {
          ori = 1;
        }
      }
      else
      {
        if ( vm < vmN )
        {
          ori = 1;
        }
        else
        {
          ori = -1;
        }
      }
    }
    else
    {
      if ( Sixidx_plan[F[Fidx[face - 1] + mf] - 1] == Sixidx_plan[Md + ne - 1] || Sixidx_plan[Md + F[Fidx[face - 1] + mf] - 1] == Sixidx_plan[Md + ne - 1] )
      {
        if ( vm < vmN )
        {
          ori = -1;
        }
        else
        {
          ori = 1;
        }
      }
      else
      {
        if ( vm < vmN )
        {
          ori = 1;
        }
        else
        {
          ori = -1;
        }
      }
    }
  }
  else
  {
    if ( Sixidx_plan[F[Fidx[face - 1] + mf] - 1] == Sixidx_plan[ne - 1] || Sixidx_plan[Md + F[Fidx[face - 1] + mf] - 1] == Sixidx_plan[ne - 1] )
    {
      if ( vm < vmN )
      {
        ori = -1;
      }
      else
      {
        ori = 1;
      }
    }
    else
    {
      if ( vm < vmN )
      {
        ori = 1;
      }
      else
      {
        ori = -1;
      }
    }
  }

  return ori;
}

int Traverse_Paths_Tree( int M, int* eparent, int* vparent, int* comb_ancestor_sparse, int* comb_ancestor_first, int* comb_ancestor_next, int z, int d, int* fin, int* pidx, int* pathedge, int* fi, int* order )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Traverse_Paths_Tree                                                    */
  /*                                                                         */
  /*  Used to process the tree created in Paths_Tree.                        */
  /*                                                                         */
  /***************************************************************************/

  int t;
  t = 1;
  int* copy_comb_ancestor_first = (int*)malloc( M * sizeof( int ) );
  int i;
  memcpy( copy_comb_ancestor_first, comb_ancestor_first, M * sizeof( int ) );
  int next;
  next = z;
  int c;
  int prev;
  int f1;
  int donep;
  donep = 0;
  int j;
  int ordcnt = 0;

  while ( t > 0 )
  {
    if ( copy_comb_ancestor_first[next - 1] != 0 )
    {
      c = comb_ancestor_sparse[copy_comb_ancestor_first[next - 1] - 1];

      if ( c > 0 )
      {
        prev = next;
        next = c;
        copy_comb_ancestor_first[prev - 1] = comb_ancestor_next[copy_comb_ancestor_first[prev - 1] - 1];
      }
      else
      {
        for ( i = 0; i < d; i++ )
        { // need d as input
          if ( fin[i] == -c )
          {
            ordcnt = ordcnt + 1;
            order[ordcnt - 1] = fi[i];
            break;
          }
        }
        prev = next;
        next = vparent[-c - 1];
        copy_comb_ancestor_first[prev - 1] = comb_ancestor_next[copy_comb_ancestor_first[prev - 1] - 1];
      }
    }
    else
    {
      if ( next == z )
      {
        break;
      }
      else
      {
        if ( comb_ancestor_first[next - 1] == 0 )
        {
          for ( i = 0; i < d; i++ )
          { // need sizes of pathedge as input. check using pathedge matrix correctly
            for ( j = pidx[i] + 1; j < pidx[i + 1] + 1; j++ )
            {
              if ( pathedge[j - 1] == next )
              {
                ordcnt = ordcnt + 1;
                order[ordcnt - 1] = fi[i];
                donep = 1;
                break;
              }
            }
            if ( donep == 1 )
            {
              donep = 0;
              break;
            }
          }
        }
        next = eparent[next - 1];
      }
    }
  }

  free( copy_comb_ancestor_first );

  return ordcnt;
}

int Paths_Tree( int M, int* pidx, int* pathedge, int* paths, int* Sixidx, int d, int* fin, int* F, int* Fdir, int* Fidx, int* pdist, int maxpdist, int* pathsidx, int* comb_ancestor_sparse, int* comb_ancestor_first, int* comb_ancestor_next, int* eparent, int* vparent, int* comb_count, int* comb_ancestor_recent )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Paths_Tree                                                             */
  /*                                                                         */
  /*  Once the destination face is found for a vertex, the paths to each of  */
  /*  its edges to its neighbours is also found. We can be certain that none */
  /*  of these edges need to cross each other, since in order to cross each  */
  /*  other they would need to reach the same face, and they can reach it on */
  /*  the same (optimal) path. However, we if two or more edges cross into   */
  /*  the same face, we need to determine which order they enter it to make  */
  /*  certain they don't end up crossing. Hence, a tree of the various paths */
  /*  is created, which can then be used cleverly in Traverse_Paths_Tree.    */
  /*                                                                         */
  /***************************************************************************/

  int i;

  for ( i = 0; i < M; i++ )
  {
    comb_ancestor_sparse[i] = 0;
    comb_ancestor_first[i] = 0;
    comb_ancestor_next[i] = 0;
    comb_count[i] = 0;
    comb_ancestor_recent[i] = 0;
  }

  int comb_ancestor_count = 0;
  int* chk = (int*)malloc( d * sizeof( int ) );
  int* checked = (int*)malloc( d * sizeof( int ) );
  for ( i = 0; i < d; i++ )
  {
    chk[i] = 0;
    checked[i] = 0;
  }
  int* stack = (int*)malloc( ( d * maxpdist + 1 ) * sizeof( int ) );
  int* edgeentry = (int*)malloc( ( d * maxpdist + 1 ) * sizeof( int ) );
  int* diststack = (int*)malloc( ( d * maxpdist + 1 ) * sizeof( int ) );
  for ( i = 0; i < d * maxpdist + 1; i++ )
  {
    stack[i] = 0;
    edgeentry[i] = 0;
    diststack[i] = 0;
  }
  int* q = (int*)malloc( d * sizeof( int ) );
  int j;
  int k;
  int chksum = 0;
  int stacksize;

  int stackpointer;
  int pn;
  int currdist;
  int en;
  int currpos;
  int f1[d];
  int fcnt = 0;
  int c;
  int c2 = 0;
  int c3 = 0;

  for ( i = 0; i < d; i++ )
  {
    if ( pdist[i] > 0 && chk[i] == 0 )
    {
      chksum = 0;
      for ( j = 0; j < d; j++ )
      {
        if ( pdist[j] > 0 && pathedge[pidx[j]] == pathedge[pidx[i]] )
        {
          chk[j] = 1;
          chksum = chksum + 1;
        }
      }
      if ( chksum > 1 )
      {
        for ( j = 0; j < d; j++ )
        {

          if ( pdist[j] > 0 && pathedge[pidx[j]] == pathedge[pidx[i]] )
          {
            stack[0] = paths[pathsidx[j] + 1];
            break;
          }
        }
        edgeentry[0] = pathedge[pidx[i]];
        stacksize = 1;
        diststack[0] = 2;
        stackpointer = 1;

        while ( stackpointer <= stacksize )
        {
          pn = stack[stackpointer - 1];
          en = edgeentry[stackpointer - 1];
          currdist = diststack[stackpointer - 1];

          stackpointer = stackpointer + 1;
          fcnt = 0;

          for ( j = Fidx[pn - 1]; j < Fidx[pn]; j++ )
          {

            if ( F[j] == en )
            {
              currpos = j;
              break;
            }
          }
          for ( j = 0; j < Fidx[pn] - Fidx[pn - 1]; j++ )
          {
            if ( currdist <= maxpdist )
            {
              fcnt = 0;
              for ( k = 0; k < d; k++ )
              {
                f1[k] = 0;
              }
              for ( k = 0; k < d; k++ )
              {
                if ( pdist[k] > 0 && pidx[k] + currdist <= pidx[k + 1] && pathedge[pidx[k] + currdist - 1] == F[currpos] )
                {
                  f1[k] = 1;
                  fcnt = fcnt + 1;
                }
              }
              if ( fcnt > 0 )
              {
                comb_count[en - 1] = comb_count[en - 1] + 1;
                comb_ancestor_count = comb_ancestor_count + 1;
                c2 = c2 + 1;
                comb_ancestor_sparse[comb_ancestor_count - 1] = F[currpos];
                if ( comb_ancestor_first[en - 1] == 0 )
                {
                  comb_ancestor_first[en - 1] = comb_ancestor_count;
                }
                if ( comb_ancestor_recent[en - 1] != 0 )
                {
                  comb_ancestor_next[comb_ancestor_recent[en - 1] - 1] = comb_ancestor_count;
                }
                comb_ancestor_recent[en - 1] = comb_ancestor_count;
                eparent[F[currpos] - 1] = en;

                if ( fcnt > 1 )
                {
                  stacksize = stacksize + 1;

                  for ( k = 0; k < d; k++ )
                  {
                    if ( f1[k] == 1 )
                    {
                      stack[stacksize - 1] = paths[pathsidx[k] + currdist];
                      break;
                    }
                  }

                  diststack[stacksize - 1] = currdist + 1;
                  edgeentry[stacksize - 1] = F[currpos];
                }
              }
            }

            c = 0;
            if ( Fdir[currpos] == 1 )
            {
              for ( k = 0; k < d; k++ )
              {
                if ( fin[k] == Sixidx[M - 1 + F[currpos]] )
                {
                  c = k + 1;
                  break;
                }
              }
            }
            else
            {
              for ( k = 0; k < d; k++ )
              {
                if ( fin[k] == Sixidx[F[currpos] - 1] )
                {
                  c = k + 1;
                  break;
                }
              }
            }
            if ( c > 0 )
            {
              if ( pdist[c - 1] == currdist - 1 && paths[pathsidx[c] - 1] == pn && pathedge[pidx[c] - 1] == en && checked[c - 1] == 0 )
              {
                checked[c - 1] = 1;
                comb_count[en - 1] = comb_count[en - 1] + 1;
                comb_ancestor_count = comb_ancestor_count + 1;
                c3 = c3 + 1;
                comb_ancestor_sparse[comb_ancestor_count - 1] = -fin[c - 1];
                if ( comb_ancestor_first[en - 1] == 0 )
                {
                  comb_ancestor_first[en - 1] = comb_ancestor_count;
                }
                if ( comb_ancestor_recent[en - 1] != 0 )
                {
                  comb_ancestor_next[comb_ancestor_recent[en - 1] - 1] = comb_ancestor_count;
                }
                comb_ancestor_recent[en - 1] = comb_ancestor_count;
                vparent[fin[c - 1] - 1] = en;
              }
            }
            if ( currpos == Fidx[pn] - 1 )
            {
              currpos = Fidx[pn - 1];
            }
            else
            {
              currpos = currpos + 1;
            }
          }
        }
      }
    }
  }

  free( chk );
  free( checked );
  free( stack );
  free( edgeentry );
  free( diststack );
  free( q );

  return comb_ancestor_count;
}

void New_Clabel( int M, int Md, int* Sixidx_plan, int* Sixidx, int* Clabel, int* Dlabel, int* Cidx, int* elist, int* eidx, int* pathedge, int* paths, int* pidx, int* pathsidx, int* pdist, int maxpdist, int d, int* fin, int* fi, int* F, int* Fdir, int* Fidx, int vm, int* origedgelabel, int* origedgeorder, int newface, int* pathedgeonce, int treecheck, int maxdeadjm )
{

  /***************************************************************************/
  /*                                                                         */
  /*  New_Clabel                                                             */
  /*                                                                         */
  /*  Once a vertex has been selected, and its new destination face chosen,  */
  /*  we need to update the data, specifically the crossing label, and also  */
  /*  Sixidx. This is done by considering the paths of each edge incident to */
  /*  the vertex which has been moved. If necessary, Paths_Tree is called to */
  /*  build a tree of movements, and Traverse_Paths_Tree is called to make   */
  /*  use of it. This only occurs if multiple edges travel to the same face. */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int change;
  int comb_ancestor_count = 0;
  int* comb_ancestor_sparse;
  comb_ancestor_sparse = (int*)malloc( Md * sizeof( int ) );
  int* comb_ancestor_first;
  comb_ancestor_first = (int*)malloc( Md * sizeof( int ) );
  int* comb_ancestor_next;
  comb_ancestor_next = (int*)malloc( Md * sizeof( int ) );
  int* eparent;
  eparent = (int*)malloc( Md * sizeof( int ) );
  int* vparent;
  vparent = (int*)malloc( Md * sizeof( int ) );
  int* comb_count;
  comb_count = (int*)malloc( Md * sizeof( int ) );
  int* comb_ancestor_recent;
  comb_ancestor_recent = (int*)malloc( Md * sizeof( int ) );
  if ( treecheck == 1 )
  {
    for ( i = 0; i < Md; i++ )
    {
      comb_ancestor_sparse[i] = 0;
      comb_ancestor_first[i] = 0;
      comb_ancestor_next[i] = 0;
      eparent[i] = 0;
      vparent[i] = 0;
      comb_count[i] = 0;
      comb_ancestor_recent[i] = 0;
    }
    comb_ancestor_count = Paths_Tree( Md, pidx, pathedge, paths, Sixidx_plan, d, fin, F, Fdir, Fidx, pdist, maxpdist, pathsidx, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, eparent, vparent, comb_count, comb_ancestor_recent );
  }

  int* newClabel;
  int* newDlabel;
  newClabel = (int*)malloc( Cidx[M] * sizeof( int ) );
  newDlabel = (int*)malloc( Cidx[M] * sizeof( int ) );
  for ( i = 0; i < Cidx[M]; i++ )
  {
    newClabel[i] = 0;
    newDlabel[i] = 0;
  }
  for ( i = 0; i < d; i++ )
  {
    change = Cidx[fi[i]] - Cidx[fi[i] - 1];
    for ( j = Cidx[fi[i]]; j < Cidx[M]; j++ )
    {
      Clabel[j - change] = Clabel[j];
      Dlabel[j - change] = Dlabel[j];
    }
    for ( j = fi[i]; j < M + 1; j++ )
    {
      Cidx[j] = Cidx[j] - change;
    }
  }

  int* newCidx = (int*)malloc( ( M + 1 ) * sizeof( int ) );
  for ( i = 0; i < M + 1; i++ )
  {
    newCidx[i] = 0;
  }
  int ci;
  int q;
  int qn;
  int k;
  int r;
  int x1cnt = -1;
  int y1cnt = -1;
  int* x1;
  x1 = (int*)malloc( d * sizeof( int ) );
  int* y1;
  y1 = (int*)malloc( d * sizeof( int ) );
  int same;
  int backe;
  int epos;
  int f1;
  int last_i = 0;
  int ordcnt;
  int* order = (int*)malloc( d * sizeof( int ) );
  int face;
  int* ctemp;
  ctemp = (int*)malloc( M * sizeof( int ) ); // max crossings size, but M for now
  for ( i = 0; i < M; i++ )
    ctemp[i] = 0;
  int ctempcnt = -1;
  int ori;

  int dtempcnt = -1;
  int* dtemp;
  dtemp = (int*)malloc( M * sizeof( int ) );

  int* fi_check;
  fi_check = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
  {
    fi_check[i] = 0;
  }
  for ( i = 0; i < d; i++ )
  {
    fi_check[fi[i] - 1] = 1;
  }

  int newccnt = -1;
  int ec;
  int curre;

  for ( i = 0; i < M; i++ )
  {
    ci = 1;
    for ( j = 0; j < d; j++ )
    {
      if ( i + 1 == fi[j] )
      {
        ci = 0;
        break;
      }
    }

    if ( ci == 1 )
    {
      ctempcnt = -1;
      dtempcnt = -1;
      x1cnt = -1;
      y1cnt = -1;
      qn = eidx[i + 1] - eidx[i];
      for ( j = 0; j < qn; j++ )
      {
        q = elist[eidx[i] + j];
        if ( treecheck == 1 && comb_count[q - 1] )
        {
          ordcnt = Traverse_Paths_Tree( Md, eparent, vparent, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, q, d, fin, pidx, pathedge, fi, order );
          x1cnt = -1;
          y1cnt = -1;
          for ( k = 0; k < d; k++ )
          {
            for ( r = pidx[k]; r < pidx[k + 1]; r++ )
            {
              if ( pathedge[r] == q )
              {
                x1cnt = x1cnt + 1;
                y1cnt = y1cnt + 1;
                x1[x1cnt] = k + 1;
                y1[y1cnt] = r - pidx[k] + 1;
              }
            }
          }

          same = 1;
          if ( origedgeorder[q - 1] == 1 || origedgeorder[q - 1] == ( eidx[origedgelabel[q - 1]] - eidx[origedgelabel[q - 1] - 1] ) )
          {
            if ( Sixidx[origedgelabel[q - 1] - 1] != Sixidx_plan[q - 1] && Sixidx[M + origedgelabel[q - 1] - 1] != Sixidx_plan[Md + q - 1] )
            {
              same = 0;
            }
          }
          else
          {
            backe = elist[eidx[origedgelabel[q - 1] - 1] + origedgeorder[q - 1] - 2];
            if ( Sixidx_plan[backe - 1] != Sixidx_plan[q - 1] && Sixidx_plan[Md + backe - 1] != Sixidx_plan[q - 1] )
            {
              same = 0;
            }
          }
          for ( k = Fidx[paths[pathsidx[x1[0] - 1] + y1[0]] - 1]; k < Fidx[paths[pathsidx[x1[0] - 1] + y1[0]]]; k++ )
          {
            if ( F[k] == q )
            {
              epos = k;
              break;
            }
          }
          if ( Fdir[epos] == 1 )
          {
            if ( same == 1 )
            {

              same = 0;
            }
            else
            {
              same = 1;
            }
          }
          if ( same == 0 )
          {
            for ( k = 0; k < ordcnt; k++ )
            {
              ctempcnt = ctempcnt + 1;
              ctemp[ctempcnt] = order[ordcnt - k - 1];
            }
            for ( k = 0; k < ordcnt; k++ )
            {
              for ( r = 0; r < d; r++ )
              {
                if ( fi[r] == order[ordcnt - k - 1] )
                {
                  f1 = r + 1;
                  break;
                }
              }
              face = paths[pathsidx[x1[x1cnt - k] - 1] + y1[y1cnt - k]];
              ori = New_Crossing_Orientation( Md, Sixidx_plan, Fidx[face] - Fidx[face - 1], F, Fidx, origedgeorder, q, vm, fin[f1 - 1], face, origedgelabel[q - 1], elist, eidx );

              dtempcnt = dtempcnt + 1;
              dtemp[dtempcnt] = ori;
            }
          }
          else
          {
            for ( k = 0; k < ordcnt; k++ )
            {
              ctempcnt = ctempcnt + 1;

              ctemp[ctempcnt] = order[k];
            }
            for ( k = 0; k < ordcnt; k++ )
            {
              for ( r = 0; r < d; r++ )
              {
                if ( fi[r] == order[k] )
                {
                  f1 = r + 1;
                  break;
                }
              }
              face = paths[pathsidx[x1[k] - 1] + y1[k]];
              ori = New_Crossing_Orientation( Md, Sixidx_plan, Fidx[face] - Fidx[face - 1], F, Fidx, origedgeorder, q, vm, fin[f1 - 1], face, origedgelabel[q - 1], elist, eidx );
              dtempcnt = dtempcnt + 1;
              dtemp[dtempcnt] = ori;
            }
          }
        }
        else
        {
          x1[0] = pathedgeonce[q - 1];
          y1[0] = pathedgeonce[maxdeadjm + q - 1]; // check this length of pathedgeonce.
          if ( x1[0] != 0 )
          {
            ctempcnt = ctempcnt + 1;
            ctemp[ctempcnt] = fi[x1[0] - 1];

            face = paths[pathsidx[x1[0] - 1] + y1[0]];
            ori = New_Crossing_Orientation( Md, Sixidx_plan, Fidx[face] - Fidx[face - 1], F, Fidx, origedgeorder, q, vm, fin[x1[0] - 1], face, origedgelabel[q - 1], elist, eidx );
            dtempcnt = dtempcnt + 1;
            dtemp[dtempcnt] = ori;
          }
        }
        if ( j < qn - 1 && fi_check[Clabel[Cidx[i] + j] - 1] == 0 )
        {
          ctempcnt = ctempcnt + 1;
          ctemp[ctempcnt] = Clabel[Cidx[i] + j];
          dtempcnt = dtempcnt + 1;
          dtemp[dtempcnt] = Dlabel[Cidx[i] + j];
        }
      }
      if ( ctempcnt != -1 )
      {
        for ( k = 0; k < ctempcnt + 1; k++ )
        {
          newccnt = newccnt + 1;
          newClabel[newccnt] = ctemp[k];
          newDlabel[newccnt] = dtemp[k];
        }
      }

      if ( last_i != 0 )
      {
        for ( k = last_i + 2; k <= i + 1; k++ )
        {
          newCidx[k - 1] = newCidx[last_i];
        }
      }
      last_i = i + 1;
      newCidx[i + 1] = newCidx[i] + ctempcnt + 1;

      if ( i != M - 1 )
      {
        for ( k = i + 3; k <= M + 1; k++ )
        {
          newCidx[k - 1] = newCidx[i + 1];
        }
      }
    }
  }
  free( ctemp );
  free( dtemp );
  free( fi_check );

  for ( i = 0; i < newCidx[M]; i++ )
  {
    Clabel[i] = newClabel[i];
    Dlabel[i] = newDlabel[i];
  }

  free( newClabel );
  free( newDlabel );

  for ( i = 0; i < M + 1; i++ )
  {
    Cidx[i] = newCidx[i];
  }

  free( newCidx );

  int nf = Fidx[newface] - Fidx[newface - 1];
  int* fincheck = (int*)malloc( d * sizeof( int ) );
  for ( i = 0; i < d; i++ )
  {
    fincheck[i] = 0;
  }
  int* atemp = (int*)malloc( d * sizeof( int ) );
  for ( i = 0; i < d; i++ )
  {
    if ( Sixidx[fi[i] - 1] == vm )
    {
      ec = Sixidx[( 4 * M ) + fi[i] - 1];
      if ( Sixidx[ec - 1] == fin[i] )
      {
        Sixidx[( 3 * M ) + ec - 1] = Sixidx[( 5 * M ) + fi[i] - 1];
      }
      else
      {
        Sixidx[( 5 * M ) + ec - 1] = Sixidx[( 5 * M ) + fi[i] - 1];
      }
      ec = Sixidx[( 5 * M ) + fi[i] - 1];
      if ( Sixidx[ec - 1] == fin[i] )
      {
        Sixidx[( 2 * M ) + ec - 1] = Sixidx[( 4 * M ) + fi[i] - 1];
      }
      else
      {
        Sixidx[( 4 * M ) + ec - 1] = Sixidx[( 4 * M ) + fi[i] - 1];
      }
    }
    else
    {
      ec = Sixidx[( 2 * M ) + fi[i] - 1];
      if ( Sixidx[ec - 1] == fin[i] )
      {
        Sixidx[( 3 * M ) + ec - 1] = Sixidx[( 3 * M ) + fi[i] - 1];
      }
      else
      {
        Sixidx[( 5 * M ) + ec - 1] = Sixidx[( 3 * M ) + fi[i] - 1];
      }
      ec = Sixidx[( 3 * M ) + fi[i] - 1];
      if ( Sixidx[ec - 1] == fin[i] )
      {
        Sixidx[( 2 * M ) + ec - 1] = Sixidx[( 2 * M ) + fi[i] - 1];
      }
      else
      {
        Sixidx[( 4 * M ) + ec - 1] = Sixidx[( 2 * M ) + fi[i] - 1];
      }
    }
  }
  int vstart;
  int acnt = -1;
  int nexteidx;
  x1cnt = -1;
  for ( i = 0; i < nf; i++ )
  {
    q = -1;
    if ( Fdir[Fidx[newface - 1] + i] == 1 )
    {
      vstart = Sixidx_plan[F[Fidx[newface - 1] + i] - 1];
    }
    else
    {
      vstart = Sixidx_plan[Md + F[Fidx[newface - 1] + i] - 1];
    }
    for ( j = 0; j < d; j++ )
    {
      if ( fin[j] == vstart )
      {
        q = j;
        break;
      }
    }

    if ( q != -1 && fincheck[q] == 0 )
    {

      acnt = acnt + 1;
      atemp[acnt] = fi[q];
      fincheck[q] = 1;
      curre = fi[q];
      if ( i > 0 )
      {
        nexteidx = Fidx[newface - 1] + i - 1;
      }
      else
      {
        nexteidx = Fidx[newface] - 1;
      }
      if ( Sixidx[origedgelabel[curre - 1] - 1] == fin[q] )
      {
        Sixidx[( 2 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[nexteidx] - 1];
        Sixidx[( 3 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[Fidx[newface - 1] + i] - 1];
      }
      else
      {
        Sixidx[( 4 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[nexteidx] - 1];
        Sixidx[( 5 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[Fidx[newface - 1] + i] - 1];
      }
      if ( Sixidx[origedgelabel[F[Fidx[newface - 1] + i] - 1] - 1] == fin[q] )
      {
        Sixidx[( 2 * M ) + origedgelabel[F[Fidx[newface - 1] + i] - 1] - 1] = fi[q];
      }
      else
      {
        Sixidx[( 4 * M ) + origedgelabel[F[Fidx[newface - 1] + i] - 1] - 1] = fi[q];
      }
      if ( Sixidx[origedgelabel[F[nexteidx] - 1] - 1] == fin[q] )
      {
        Sixidx[( 3 * M ) + origedgelabel[F[nexteidx] - 1] - 1] = fi[q];
      }
      else
      {
        Sixidx[( 5 * M ) + origedgelabel[F[nexteidx] - 1] - 1] = fi[q];
      }
    }
    x1cnt = -1;
    for ( j = 0; j < d; j++ )
    {
      if ( pdist[j] > 0 && pathedge[pidx[j]] == F[Fidx[newface - 1] + i] )
      {
        x1cnt = x1cnt + 1;
        x1[x1cnt] = j + 1;
      }
    }

    if ( x1cnt > 0 )
    {
      ordcnt = Traverse_Paths_Tree( Md, eparent, vparent, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, F[Fidx[newface - 1] + i], d, fin, pidx, pathedge, fi, order );
      for ( j = 0; j < ordcnt; j++ )
      {
        acnt = acnt + 1;
        atemp[acnt] = order[j];
      }
    }
    else if ( x1cnt == 0 )
    {
      acnt = acnt + 1;
      atemp[acnt] = fi[x1[0] - 1];
    }
  }

  free( fincheck );

  int fs;
  int bs;
  for ( i = 0; i < d; i++ )
  {
    if ( i == 0 )
    {
      bs = atemp[acnt];
    }
    else
    {
      bs = atemp[i - 1];
    }
    if ( i == d - 1 )
    {
      fs = atemp[0];
    }
    else
    {
      fs = atemp[i + 1];
    }

    if ( Sixidx[atemp[i] - 1] == vm )
    {
      Sixidx[( 2 * M ) + atemp[i] - 1] = fs;
      Sixidx[( 3 * M ) + atemp[i] - 1] = bs;
    }
    else
    {
      Sixidx[( 4 * M ) + atemp[i] - 1] = fs;
      Sixidx[( 5 * M ) + atemp[i] - 1] = bs;
    }
  }
  free( atemp );

  int* x12 = (int*)malloc( maxpdist * sizeof( int ) );
  int* checked;
  int* newedgelabel;
  int* newedgelabelun;
  if ( maxpdist > 0 )
  {
    checked = (int*)malloc( maxpdist * sizeof( int ) );
    newedgelabel = (int*)malloc( maxpdist * sizeof( int ) );
    newedgelabelun = (int*)malloc( maxpdist * sizeof( int ) );
  }
  int tj;
  int sr;
  for ( i = 0; i < d; i++ )
  {
    fs = -1;
    if ( pdist[i] > 0 )
    {
      x1cnt = -1;
      fs = -1;
      for ( j = Cidx[fi[i] - 1]; j < Cidx[M]; j++ )
      {
        fs = fs + 1;
        Clabel[Cidx[M] - 1 + pdist[i] - fs] = Clabel[Cidx[M] - 1 - fs];
        Dlabel[Cidx[M] - 1 + pdist[i] - fs] = Dlabel[Cidx[M] - 1 - fs];
      }
      fs = -1;
      if ( vm < fin[i] )
      {
        for ( j = Cidx[fi[i] - 1]; j < Cidx[fi[i] - 1] + pdist[i]; j++ )
        {
          fs = fs + 1;
          Clabel[j] = origedgelabel[pathedge[pidx[i] + fs] - 1];
        }
      }
      else
      {
        for ( j = Cidx[fi[i] - 1]; j < Cidx[fi[i] - 1] + pdist[i]; j++ )
        {
          fs = fs + 1;
          Clabel[j] = origedgelabel[pathedge[pidx[i + 1] - 1 - fs] - 1];
        }
      }
      for ( j = fi[i]; j <= M; j++ )
      {
        Cidx[j] = Cidx[j] + pdist[i];
      }
      for ( j = 0; j < maxpdist; j++ )
      {
        checked[j] = 0;
      }
      for ( j = 0; j < pdist[i]; j++ )
      {
        x1cnt = -1;
        for ( k = 0; k < pdist[i]; k++ )
        {
          if ( Clabel[Cidx[fi[i] - 1] + j] == Clabel[Cidx[fi[i] - 1] + k] )
          {
            x1cnt = x1cnt + 1;
            x12[x1cnt] = k + 1;
          }
        }
        if ( x1cnt > 0 && checked[j] == 0 )
        {

          for ( k = 0; k <= x1cnt; k++ )
          {
            checked[x12[k] - 1] = 1;
            if ( vm < fin[i] )
            {
              newedgelabel[k] = pathedge[pidx[i] + x12[k] - 1];
            }
            else
            {
              newedgelabel[k] = pathedge[pidx[i + 1] - x12[k]];
            }
            newedgelabelun[k] = newedgelabel[k];
          }
          qsort( newedgelabel, x1cnt + 1, sizeof( int ), ascend_cmpfunc );
          tj = Clabel[Cidx[fi[i] - 1] + j];
          fs = -1;
          for ( k = Cidx[tj - 1]; k < Cidx[tj]; k++ )
          {
            if ( Clabel[k] == fi[i] )
            {
              fs = fs + 1;
              sr = -1;
              for ( r = 0; r < x1cnt + 1; r++ )
              {
                if ( newedgelabelun[r] == newedgelabel[fs] )
                {
                  sr = r;
                }
              }
              Dlabel[Cidx[fi[i] - 1] + x12[sr] - 1] = -Dlabel[k]; // may be wrong here.
            }
          }
        }
        else if ( x1cnt == 0 )
        {
          tj = Clabel[Cidx[fi[i] - 1] + j];
          for ( k = Cidx[tj - 1]; k < Cidx[tj]; k++ )
          {

            if ( Clabel[k] == fi[i] )
            {
              Dlabel[Cidx[fi[i] - 1] + j] = -Dlabel[k];
              break;
            }
          }
        }
      }
    }
  }
  free( x12 );
  free( x1 );
  free( y1 );
  if ( maxpdist > 0 )
  {
    free( newedgelabel );
    free( newedgelabelun );
    free( checked );
  }

  int* checked2 = (int*)malloc( d * sizeof( int ) );
  for ( j = 0; j < d; j++ )
  {
    checked2[j] = 0;
  }
  int* endfaces = (int*)malloc( d * sizeof( int ) );
  int startedge;
  int starte;
  int fn;
  int count;
  int fk;
  int ccnt = -1;
  int* ctemp2 = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
    ctemp2[i] = 150;
  int facee;

  for ( j = 0; j < d; j++ )
    endfaces[j] = paths[pathsidx[j + 1] - 1];

  for ( i = 0; i < d; i++ )
  {
    if ( checked2[i] == 0 )
    {
      for ( j = 0; j < d; j++ )
      {
        if ( endfaces[j] == endfaces[i] )
        {
          checked2[j] = 1;
        }
      }
      if ( pdist[i] > 0 )
      {
        startedge = pathedge[pidx[i + 1] - 1];
        facee = endfaces[i];
        ccnt = -1;

        if ( treecheck > 0 && comb_count[startedge - 1] > 0 )
        {
          ordcnt = Traverse_Paths_Tree( Md, eparent, vparent, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, startedge, d, fin, pidx, pathedge, fi, order );
          for ( k = 0; k < ordcnt; k++ )
          {
            ccnt = ccnt + 1;
            ctemp2[ccnt] = order[k];
          }
        }
        else
        {
          ccnt = ccnt + 1;
          ctemp2[ccnt] = fi[i];
        }
        for ( j = Fidx[facee - 1]; j < Fidx[facee]; j++ )
        {
          if ( F[j] == startedge )
          {
            fs = j;
            break;
          }
        }
        fn = Fidx[facee] - Fidx[facee - 1];
        count = 1;
        for ( j = 0; j < fn; j++ )
        {
          if ( fs < Fidx[facee] - 1 )
          {
            fs = fs + 1;
          }
          else
          {
            fs = Fidx[facee - 1];
          }

          if ( fs > Fidx[facee - 1] )
          {
            bs = fs - 1;
          }
          else
          {
            bs = Fidx[facee] - 1;
          }

          if ( Fdir[fs] == 1 )
          {
            vstart = Sixidx_plan[F[fs] - 1];
          }
          else
          {
            vstart = Sixidx_plan[Md + F[fs] - 1];
          }
          x1cnt = 0;
          fk = -1;
          for ( k = 0; k < d; k++ )
          {
            if ( fin[k] == vstart )
            {
              fk = k;
              x1cnt = x1cnt + 1;
              if ( endfaces[k] != facee )
              {
                x1cnt = -1;
                break;
              }
            }
          }
          if ( fk != -1 && x1cnt > 0 )
          {
            if ( ctemp2[count - 1] == fi[fk] )
            {
              count = count + 1;
              if ( Sixidx[fi[fk] - 1] == vstart )
              {
                Sixidx[( 2 * M ) + fi[fk] - 1] = origedgelabel[F[bs] - 1];
                Sixidx[( 3 * M ) + fi[fk] - 1] = origedgelabel[F[fs] - 1];
              }
              else
              {
                Sixidx[( 4 * M ) + fi[fk] - 1] = origedgelabel[F[bs] - 1];
                Sixidx[( 5 * M ) + fi[fk] - 1] = origedgelabel[F[fs] - 1];
              }
              if ( Sixidx[origedgelabel[F[fs] - 1] - 1] == vstart )
              {
                Sixidx[( 2 * M ) + origedgelabel[F[fs] - 1] - 1] = fi[fk];
              }
              else
              {
                Sixidx[( 4 * M ) + origedgelabel[F[fs] - 1] - 1] = fi[fk];
              }
              if ( Sixidx[origedgelabel[F[bs] - 1] - 1] == vstart )
              {
                Sixidx[( 3 * M ) + origedgelabel[F[bs] - 1] - 1] = fi[fk];
              }
              else
              {
                Sixidx[( 5 * M ) + origedgelabel[F[bs] - 1] - 1] = fi[fk];
              }
            }
          }
          x1cnt = 0;
          for ( k = 0; k < pidx[d]; k++ )
          {
            if ( pathedge[k] == F[fs] )
            {
              x1cnt = x1cnt + 1;
            }
          }
          if ( x1cnt > 0 && F[fs] != startedge )
          {
            count = count + x1cnt;
          }
          if ( count > ccnt + 1 )
          {

            break;
          }
        }
      }
    }
  }

  free( endfaces );
  free( ctemp2 );
  free( comb_ancestor_sparse );
  free( comb_ancestor_first );
  free( comb_ancestor_next );
  free( eparent );
  free( vparent );
  free( comb_count );
  free( checked2 );
  free( comb_ancestor_recent );
  free( order );
}

void Find_Merge_Faces( int Nd, int M_plan, int* F, int* Fdir, int* Fidx, int* dualel, bool* bad_edges, int* paths, int* pathsidx, int* pdist, int d, int oldNd, int* pooledgemap, int* newf, int* newfidx, int* newdir, int* rf_map, int maxmaxpaths )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Find_Merge_Faces                                                       */
  /*                                                                         */
  /*  Once an improvement has been found, we need to map the merged faces    */
  /*  back to the original faces, so that New_Clabel can be executed.        */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;
  int cnt3 = 0;
  int newfnum = 0;

  for ( i = 0; i < Fidx[Nd]; i++ )
  {
    newf[i] = 0;
    newdir[i] = 0;
  }

  for ( i = 0; i < Nd; i++ )
    newfidx[i] = 0;

  for ( i = 0; i < maxmaxpaths; i++ )
    rf_map[i] = 0;

  for ( i = 1; i <= d; i++ )
  {
    for ( j = 1; j <= pdist[i - 1] + 1; j++ )
    {
      if ( rf_map[paths[pathsidx[i - 1] + j - 1] - 1] == 0 )
        if ( paths[pathsidx[i - 1] + j - 1] <= oldNd )
        {
          newfnum++;
          newfidx[newfnum - 1] = cnt3;
          int f = paths[pathsidx[i - 1] + j - 1];
          rf_map[f - 1] = newfnum;
          for ( k = 1; k <= Fidx[f] - Fidx[f - 1]; k++ )
          {
            newf[cnt3 + k - 1] = F[Fidx[f - 1] + k - 1];
            newdir[cnt3 + k - 1] = Fdir[Fidx[f - 1] + k - 1];
          }
          cnt3 = cnt3 + Fidx[f] - Fidx[f - 1];
        }
        else
        {
          newfnum++;
          rf_map[paths[pathsidx[i - 1] + j - 1] - 1] = newfnum;
          int cnt2;
          int currf = pooledgemap[paths[pathsidx[i - 1] + j - 1] - oldNd - 1];

          int starte;
          int startd;

          int stackmax = 10;
          int* stack;
          stack = (int*)malloc( stackmax * sizeof( int ) );
          stack[0] = currf;
          int stackpos = 1;
          int stackcount = 0;
          int t = 1;

          while ( t == 1 )
          {
            stackcount++;
            currf = stack[stackcount - 1];
            for ( k = Fidx[currf - 1]; k < Fidx[currf]; k++ )
            {
              if ( bad_edges[F[k] - 1] == 1 )
              {
                stackpos++;
                if ( stackpos > stackmax )
                {
                  stack = (int*)realloc( stack, 2 * stackmax * sizeof( int ) );
                  stackmax = 2 * stackmax;
                }
                if ( dualel[F[k] - 1] == currf )
                {
                  stack[stackpos - 1] = dualel[F[k] - 1 + M_plan];
                }
                else
                {
                  stack[stackpos - 1] = dualel[F[k] - 1];
                }
              }
              else
              {
                starte = F[k];
                startd = Fdir[k];
                cnt2 = k - Fidx[currf - 1] + 1;
                t = 0;
                break;
              }
            }
          }
          free( stack );
          newfidx[newfnum - 1] = cnt3;
          cnt3++;

          newf[cnt3 - 1] = starte;
          newdir[cnt3 - 1] = startd;
          while ( true )
          {
            if ( cnt2 == Fidx[currf] - Fidx[currf - 1] )
              cnt2 = 1;
            else
              cnt2++;

            int curre = F[Fidx[currf - 1] + cnt2 - 1];
            int currd = Fdir[Fidx[currf - 1] + cnt2 - 1];

            if ( bad_edges[curre - 1] )
            {
              currd = -currd;
              if ( dualel[curre - 1] != dualel[curre - 1 + M_plan] )
              {
                if ( currf == dualel[curre - 1] )
                  currf = dualel[curre - 1 + M_plan];
                else
                  currf = dualel[curre - 1];

                for ( k = Fidx[currf - 1] + 1; k <= Fidx[currf]; k++ )
                  if ( F[k - 1] == curre )
                  {
                    cnt2 = k - Fidx[currf - 1];
                    break;
                  }
              }
              else
              {
                int cnt4_1 = -1;
                int cnt4_2 = -1;
                for ( k = Fidx[currf - 1] + 1; k <= Fidx[currf]; k++ )
                  if ( F[k - 1] == curre )
                  {
                    if ( cnt4_1 == -1 )
                      cnt4_1 = k - Fidx[currf - 1];
                    else
                    {
                      cnt4_2 = k - Fidx[currf - 1];
                      break;
                    }
                  }

                if ( Fdir[Fidx[currf - 1] + cnt4_1 - 1] != currd )
                  cnt2 = cnt4_1;
                else
                {
                  if ( Fdir[Fidx[currf - 1] + cnt4_2 - 1] != currd )
                    cnt2 = cnt4_2;
                }
              }
              continue;
            }

            if ( curre == starte && currd == startd )
            {
              newfidx[newfnum] = cnt3;
              break;
            }
            cnt3++;
            newf[cnt3 - 1] = curre;
            newdir[cnt3 - 1] = currd;
          }
        }
    }
  }

  newfidx[newfnum] = cnt3;
}

void Sp_Bigface( int* dvadj, int* Fidx, int source, int N, int* facelist, int facelistrows, int nav, int* dist, int* prev )
{

  /***************************************************************************/
  /*                                                                         */
  /*  So_Bigface                                                             */
  /*                                                                         */
  /*  This is a modified version of the BFS shortest path algorithm used in  */
  /*  Get_Distances. The primary difference is that, for the bigface scheme, */
  /*  it is often possible to terminate the algorithm early, once all of the */
  /*  neighbours of the moving vertex have been reached.                     */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;

  bool* dest = (bool*)malloc( N * nav * sizeof( bool ) );
  bool* lookdest = (bool*)malloc( N * sizeof( bool ) );

  for ( i = 0; i < N * nav; i++ )
    dest[i] = false;

  for ( i = 0; i < N; i++ )
    lookdest[i] = false;

  for ( i = 0; i < facelistrows; i++ )
  {
    for ( j = 0; j < nav; j++ )
    {
      if ( facelist[i + j * facelistrows] > 0 )
      {
        lookdest[facelist[i + j * facelistrows] - 1] = true;
        dest[facelist[i + j * facelistrows] - 1 + j * N] = true;
      }
    }
  }

  bool* done = (bool*)malloc( N * sizeof( bool ) );
  for ( j = 0; j < N; j++ )
    done[j] = false;

  done[source - 1] = true;

  for ( j = 0; j < N; j++ )
    dist[j] = N;
  dist[source - 1] = 0;

  int* looknext = (int*)malloc( ( N + 1 ) * sizeof( int ) );
  for ( j = 0; j < N + 1; j++ )
    looknext[j] = 0;
  looknext[0] = source;
  int looknextcount = 1;
  int index = 0;
  int count = 0;

  int destfaces = 0;
  bool* checkfaces = (bool*)malloc( nav * sizeof( bool ) );
  for ( j = 0; j < nav; j++ )
    checkfaces[j] = false;

  while ( looknext[index] != 0 )
  {
    int v = looknext[index];

    for ( i = Fidx[v - 1]; i < Fidx[v]; i++ )
    {
      int check = dvadj[i];
      if ( done[check - 1] == false )
      {
        looknext[looknextcount] = check;
        looknextcount = looknextcount + 1;
        dist[check - 1] = dist[v - 1] + 1;
        prev[check - 1] = v;
        done[check - 1] = true;
        if ( lookdest[check - 1] == true )
        {
          for ( k = 1; k <= nav; k++ )
          {
            if ( dest[check - 1 + ( k - 1 ) * N] == true )
            {
              if ( checkfaces[k - 1] == false )
              {
                checkfaces[k - 1] = true;
                destfaces = destfaces + 1;
              }
            }
          }
          if ( destfaces == nav )
          {
            index = N - 1;
            break;
          }
        }
      }
    }
    index = index + 1;
  }

  free( dest );
  free( lookdest );
  free( done );
  free( looknext );
  free( checkfaces );
}

int Get_Distances( int* dvadjm, int* Fidxm, int Ndm, int Mdm, int M_plan, int* Sixidx, int* dualel, int requiredimprovement, int d, int* fin, int* incidente, int* deg, bool* bad_edges, int current_crossings_best, int* deadjm, int mcnt, int* merged_face_labels, int max_F_index, int scheme, int bigface_depth, int maxdegfin, int maxdeadjm, int** pathedge, int* pidx, int* pathedgeonce, int** paths, int* pathsidx, int* pdist, int* delv, int* outputs, int bigface_fails )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Get_Distances                                                          */
  /*                                                                         */
  /*  This subroutine is the biggest and most computationally expensive part */
  /*  of the algorithm. There are two sub-algorithms within, the Bigface     */
  /*  method and the Standard method.                                        */
  /*                                                                         */
  /*                                                                         */
  /*  Bigface method:   (only used if min_scheme = 1)                        */
  /*                                                                         */
  /*  In this method, the "biggest" face (that is, the face with the most    */
  /*  sides) is looked at. This will have been determined previously. This   */
  /*  face is considered as a destination for the moved vertex, and shortest */
  /*  paths for each of its neighbours are computed using Sp_Bigface. If an  */
  /*  improvement is found, it is taken. If not, the faces surrounding the   */
  /*  biggest face are considered as well, to a depth of bigface_depth. If   */
  /*  no improvement is found, the Standard method is executed.              */
  /*                                                                         */
  /*                                                                         */
  /*  Standard method:                                                       */
  /*                                                                         */
  /*  In this method, the intention is to find a destination face for the    */
  /*  moving vertex. However, finding the shortest paths from all possible   */
  /*  destination faces is too inefficient. Instead, we want to compute the  */
  /*  shortest paths from the faces surrounding each neighbour of the moving */
  /*  vertex. Any of the faces surrounding a neighbour is sufficient to get  */
  /*  to the neighbour, hence we use the following approach:                 */
  /*                                                                         */
  /*    For each neighbour of the moving vertex, consider a new graph        */
  /*    where all faces around that vertex are merged into one. In the       */
  /*    dual graph, this corresponds to contracting vertices. Then, find     */
  /*    the shortest paths (using BFS) to all other vertices in the dual     */
  /*    graph. The shortest paths to the merged vertices is taken as 0.      */
  /*                                                                         */
  /*    Then, add up all of the combined distances to the faces, and pick    */
  /*    the smallest. If it is small enough to result in fewer crossings,    */
  /*    then an improvement has been found. If not, then the moving vertex   */
  /*    is already located in its optimal face for the current iteration.    */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;

  int bigface_success = -1;

  for ( i = 0; i < 2 * maxdeadjm; i++ )
    pathedgeonce[i] = 0;

  for ( i = 0; i < d; i++ )
    pdist[i] = 0;

  for ( i = 0; i < d + 1; i++ )
  {
    pidx[i] = 0;
    pathsidx[i] = 0;
  }

  for ( i = 0; i < d * maxdegfin; i++ )
    delv[i] = 0;

  int newcross = Mdm;
  int newface = 0;
  int treecheck = 0;
  int improvementfound = 0;
  int used_bigface = 0;

  bool continue_looking = true;

  int maxpdist = 0;

  if ( scheme == 1 )
  {
    int* faces_of_interest = (int*)malloc( d * maxdegfin * sizeof( int ) );
    int* bad_face = (int*)malloc( d * maxdegfin * sizeof( int ) );
    for ( i = 0; i < d * maxdegfin; i++ )
    {
      faces_of_interest[i] = 1;
      bad_face[i] = 0;
    }

    for ( i = 1; i <= d; i++ )
    {
      for ( j = deg[fin[i - 1] - 1] + 1; j <= maxdegfin; j++ )
        bad_face[j - 1 + ( i - 1 ) * maxdegfin] = 1;

      int edge = incidente[fin[i - 1] - 1];
      for ( j = 1; j <= deg[fin[i - 1] - 1]; j++ )
      {
        if ( bad_edges[edge - 1] )
        {
          if ( Sixidx[edge - 1] == fin[i - 1] )
            edge = Sixidx[edge - 1 + 2 * M_plan];
          else
            edge = Sixidx[edge - 1 + 4 * M_plan];
          bad_face[j - 1 + ( i - 1 ) * maxdegfin] = 1;
          continue;
        }
        if ( Sixidx[edge - 1] == fin[i - 1] )
        {
          faces_of_interest[j - 1 + ( i - 1 ) * maxdegfin] = dualel[edge - 1];
          edge = Sixidx[edge - 1 + 2 * M_plan];
        }
        else
        {
          faces_of_interest[j - 1 + ( i - 1 ) * maxdegfin] = dualel[edge - 1 + M_plan];
          edge = Sixidx[edge - 1 + 4 * M_plan];
        }
      }
    }

    int* a = (int*)malloc( Ndm * sizeof( int ) );
    int* prev = (int*)malloc( Ndm * sizeof( int ) );
    for ( i = 0; i < Ndm; i++ )
    {
      a[i] = 0;
      prev[i] = 0;
    }

    int* input_faces = (int*)malloc( d * maxdegfin * sizeof( int ) );
    for ( i = 0; i < d * maxdegfin; i++ )
      input_faces[i] = faces_of_interest[i] * ( 1 - bad_face[i] );

    Sp_Bigface( dvadjm, Fidxm, max_F_index, Ndm, input_faces, maxdegfin, d, a, prev );
    for ( i = 1; i <= mcnt; i++ )
      a[merged_face_labels[i - 1] - 1] = Ndm;

    int* dist = (int*)malloc( d * sizeof( int ) );
    int* face_choices = (int*)malloc( d * sizeof( int ) );
    for ( i = 0; i < d; i++ )
    {
      dist[i] = Mdm;
      face_choices[i] = 0;
    }
    for ( i = 1; i <= d; i++ )
      for ( j = 1; j <= maxdegfin; j++ )
        if ( dist[i - 1] > a[faces_of_interest[j - 1 + ( i - 1 ) * maxdegfin] - 1] + Mdm * bad_face[j - 1 + ( i - 1 ) * maxdegfin] )
        {
          dist[i - 1] = a[faces_of_interest[j - 1 + ( i - 1 ) * maxdegfin] - 1];
          face_choices[i - 1] = j;
        }

    newcross = 0;
    for ( i = 1; i <= d; i++ )
      newcross = newcross + dist[i - 1];
    newface = max_F_index;

    if ( newcross + requiredimprovement < current_crossings_best )
    {
      continue_looking = false;
      bigface_fails = 0;
      improvementfound = 1;
      for ( i = 1; i <= d; i++ )
        pdist[i - 1] = dist[i - 1];

      for ( i = 1; i <= d; i++ )
        if ( pdist[i - 1] > maxpdist )
          maxpdist = pdist[i - 1];

      pidx[0] = 0;
      pathsidx[0] = 0;
      for ( i = 1; i <= d; i++ )
      {
        pidx[i] = pidx[i - 1] + pdist[i - 1];
        pathsidx[i] = pathsidx[i - 1] + pdist[i - 1] + 1;
      }
      *paths = (int*)malloc( pathsidx[d] * sizeof( int ) );
      for ( i = 0; i < pathsidx[d]; i++ )
        ( *paths )[i] = 0;
      if ( maxpdist > 0 )
      {
        *pathedge = (int*)malloc( pidx[d] * sizeof( int ) );
        for ( i = 0; i < pidx[d]; i++ )
          ( *pathedge )[i] = 0;
      }
      else
      {
        *pathedge = (int*)malloc( d * sizeof( int ) );
        for ( i = 0; i < d; i++ )
          ( *pathedge )[i] = 0;
      }

      bool* checked = (bool*)malloc( Ndm * sizeof( bool ) );
      for ( i = 0; i < Ndm; i++ )
        checked[i] = false;
      for ( i = 1; i <= d; i++ )
      {
        int face = faces_of_interest[face_choices[i - 1] - 1 + ( i - 1 ) * maxdegfin];
        ( *paths )[pathsidx[i - 1] + dist[i - 1]] = face;
        if ( checked[face - 1] && face != max_F_index )
          treecheck = 1;
        checked[face - 1] = true;
        for ( j = 1; j <= dist[i - 1]; j++ )
        {
          face = prev[face - 1];
          ( *paths )[pathsidx[i - 1] + dist[i - 1] - j] = face;
          if ( checked[face - 1] && face != max_F_index )
            treecheck = 1;
          checked[face - 1] = 1;
        }
        for ( j = 1; j <= dist[i - 1]; j++ )
        {
          int f1 = 0;
          for ( k = Fidxm[( *paths )[pathsidx[i - 1] + ( j - 1 )] - 1] + 1; k <= Fidxm[( *paths )[pathsidx[i - 1] + ( j - 1 )]]; k++ )
            if ( dvadjm[k - 1] == ( *paths )[pathsidx[i - 1] + j] && bad_edges[deadjm[k - 1] - 1] == false )
            {
              f1 = k;
              break;
            }

          ( *pathedge )[pidx[i - 1] + ( j - 1 )] = deadjm[f1 - 1];
          pathedgeonce[deadjm[f1 - 1] - 1] = i;
          pathedgeonce[deadjm[f1 - 1] - 1 + maxdeadjm] = j;
        }
      }
      free( checked );
    }
    free( dist );
    free( face_choices );

    if ( continue_looking && bigface_depth > 0 )
    {
      int* old_a = (int*)malloc( Ndm * sizeof( int ) );
      for ( i = 1; i <= Ndm; i++ )
        old_a[i - 1] = a[i - 1];
      int faces;
      for ( faces = 1; faces < Ndm; faces++ )
      {
        if ( faces != max_F_index && old_a[faces - 1] <= bigface_depth )
        {
          bigface_success = old_a[faces - 1];
          Sp_Bigface( dvadjm, Fidxm, faces, Ndm, input_faces, maxdegfin, d, a, prev );

          for ( i = 1; i <= mcnt; i++ )
            a[merged_face_labels[i - 1] - 1] = Mdm;

          int* dist = (int*)malloc( d * sizeof( int ) );
          int* face_choices = (int*)malloc( d * sizeof( int ) );
          for ( i = 0; i < d; i++ )
          {
            dist[i] = Mdm;
            face_choices[i] = 0;
          }

          for ( i = 1; i <= d; i++ )
            for ( j = 1; j <= maxdegfin; j++ )
              if ( dist[i - 1] > a[faces_of_interest[j - 1 + ( i - 1 ) * maxdegfin] - 1] + Mdm * bad_face[j - 1 + ( i - 1 ) * maxdegfin] )
              {
                dist[i - 1] = a[faces_of_interest[j - 1 + ( i - 1 ) * maxdegfin] - 1];
                face_choices[i - 1] = j;
              }

          newcross = 0;
          for ( i = 1; i <= d; i++ )
            newcross = newcross + dist[i - 1];
          newface = faces;
          if ( newcross + requiredimprovement < current_crossings_best )
          {
            continue_looking = false;
            bigface_fails = 0;
            improvementfound = 1;
            for ( i = 1; i <= d; i++ )
              pdist[i - 1] = dist[i - 1];

            maxpdist = 0;
            for ( i = 1; i <= d; i++ )
              if ( pdist[i - 1] > maxpdist )
                maxpdist = pdist[i - 1];

            pidx[0] = 0;
            pathsidx[0] = 0;
            for ( i = 1; i <= d; i++ )
            {
              pidx[i] = pidx[i - 1] + pdist[i - 1];
              pathsidx[i] = pathsidx[i - 1] + pdist[i - 1] + 1;
            }
            *paths = (int*)malloc( pathsidx[d] * sizeof( int ) );
            for ( i = 0; i < pathsidx[d]; i++ )
              ( *paths )[i] = 0;
            if ( maxpdist > 0 )
            {
              *pathedge = (int*)malloc( pidx[d] * sizeof( int ) );
              for ( i = 0; i < pidx[d]; i++ )
                ( *pathedge )[i] = 0;
            }
            else
            {
              *pathedge = (int*)malloc( d * sizeof( int ) );
              for ( i = 0; i < d; i++ )
                ( *pathedge )[i] = 0;
            }

            bool* checked = (bool*)malloc( Ndm * sizeof( bool ) );
            for ( i = 0; i < Ndm; i++ )
              checked[i] = false;

            for ( i = 1; i <= d; i++ )
            {
              int face = faces_of_interest[face_choices[i - 1] - 1 + ( i - 1 ) * maxdegfin];
              ( *paths )[pathsidx[i - 1] + ( dist[i - 1] )] = face;
              if ( checked[face - 1] && face != faces )
                treecheck = 1;
              checked[face - 1] = true;
              for ( j = 1; j <= dist[i - 1]; j++ )
              {
                face = prev[face - 1];
                ( *paths )[pathsidx[i - 1] + ( dist[i - 1] - j )] = face;
                if ( checked[face - 1] && face != faces )
                  treecheck = 1;
                checked[face - 1] = 1;
              }

              for ( j = 1; j <= dist[i - 1]; j++ )
              {
                int f1 = 0;
                for ( k = Fidxm[( *paths )[pathsidx[i - 1] + ( j - 1 )] - 1] + 1; k <= Fidxm[( *paths )[pathsidx[i - 1] + ( j - 1 )]]; k++ )
                  if ( dvadjm[k - 1] == ( *paths )[pathsidx[i - 1] + j] && bad_edges[deadjm[k - 1] - 1] == false )
                  {
                    f1 = k;
                    break;
                  }

                ( *pathedge )[pidx[i - 1] + ( j - 1 )] = deadjm[f1 - 1];
                pathedgeonce[deadjm[f1 - 1] - 1] = i;
                pathedgeonce[deadjm[f1 - 1] - 1 + maxdeadjm] = j;
              }
            }
            free( checked );
          }
          free( dist );
          free( face_choices );
        }
        if ( continue_looking == false )
        {
          break;
        }
      }
      free( old_a );
    }
    else
    {
      bigface_success = 0;
    }
    if ( continue_looking )
    {
      bigface_fails++;
      if ( bigface_fails >= 20 )
      {
        scheme = 2;
      }
    }

    free( bad_face );
    free( faces_of_interest );
    free( a );
    free( prev );
    free( input_faces );
  }

  if ( continue_looking )
  {
    bigface_success = -1;
    int Nd = Ndm + 1;
    int* dist = (int*)malloc( d * Nd * sizeof( int ) );
    int* prev = (int*)malloc( d * Nd * sizeof( int ) );
    for ( i = 0; i < d * Nd; i++ )
    {
      dist[i] = Nd;
      prev[i] = 0;
    }

    for ( k = 1; k <= d; k++ )
    {
      int* dvadj = (int*)malloc( 3 * Mdm * sizeof( int ) );
      int kkk;
      for ( kkk = 0; kkk < 2 * Mdm; kkk++ )
        dvadj[kkk] = dvadjm[kkk];
      for ( kkk = 2 * Mdm; kkk < 3 * Mdm; kkk++ )
        dvadj[kkk] = 0;

      int* Fidx = (int*)malloc( ( Ndm + 2 ) * sizeof( int ) );
      for ( kkk = 0; kkk < Ndm + 1; kkk++ )
        Fidx[kkk] = Fidxm[kkk];

      int curre = incidente[fin[k - 1] - 1];

      while ( bad_edges[curre - 1] )
      {
        if ( Sixidx[curre - 1] == fin[k - 1] )
          curre = Sixidx[curre - 1 + 2 * M_plan];
        else
          curre = Sixidx[curre - 1 + 4 * M_plan];
      }

      int nexte = -1;
      if ( Sixidx[curre - 1] == fin[k - 1] )
        nexte = Sixidx[curre - 1 + 2 * M_plan];
      else
        nexte = Sixidx[curre - 1 + 4 * M_plan];

      while ( bad_edges[nexte - 1] )
      {
        if ( Sixidx[nexte - 1] == fin[k - 1] )
          nexte = Sixidx[nexte - 1 + 2 * M_plan];
        else
          nexte = Sixidx[nexte - 1 + 4 * M_plan];
      }

      int* newv = (int*)malloc( Mdm * sizeof( int ) );
      for ( i = 0; i < Mdm; i++ )
        newv[i] = 0;
      int cnt = 0;

      int temp_curre = curre;
      int temp_nexte = nexte;

      for ( i = 1; i <= deg[fin[k - 1] - 1]; i++ )
      {
        int fc = -1;
        if ( Sixidx[temp_curre - 1] == fin[k - 1] )
          fc = dualel[temp_curre - 1];
        else
          fc = dualel[temp_curre - 1 + M_plan];

        delv[k - 1 + ( i - 1 ) * d] = fc;
        dist[k - 1 + ( fc - 1 ) * d] = 0;

        temp_curre = temp_nexte;

        if ( Sixidx[temp_curre - 1] == fin[k - 1] )
          temp_nexte = Sixidx[temp_curre - 1 + 2 * M_plan];
        else
          temp_nexte = Sixidx[temp_curre - 1 + 4 * M_plan];

        while ( bad_edges[temp_nexte - 1] )
        {
          if ( Sixidx[temp_nexte - 1] == fin[k - 1] )
            temp_nexte = Sixidx[temp_nexte - 1 + 2 * M_plan];
          else
            temp_nexte = Sixidx[temp_nexte - 1 + 4 * M_plan];
        }
      }
      for ( i = 1; i <= deg[fin[k - 1] - 1]; i++ )
      {

        int fc = -1;
        if ( Sixidx[curre - 1] == fin[k - 1] )
          fc = dualel[curre - 1];
        else
          fc = dualel[curre - 1 + M_plan];

        for ( j = 1; j <= Fidx[fc] - Fidx[fc - 1]; j++ )
        {

          bool createEdge = true;

          int jj;
          for ( jj = 1; jj <= maxdegfin; jj++ )
          {
            if ( delv[k - 1 + ( jj - 1 ) * d] == dvadj[Fidx[fc - 1] + j - 1] )
            {
              createEdge = false;
              break;
            }
          }

          if ( createEdge )
          {
            cnt = cnt + 1;
            newv[cnt - 1] = dvadj[Fidx[fc - 1] + j - 1];
            int* f1 = (int*)malloc( ( 8 * ( Fidx[newv[cnt - 1]] - Fidx[newv[cnt - 1] - 1] ) ) * sizeof( int ) );
            int f1cnt = 0;
            int r;
            for ( r = Fidx[newv[cnt - 1] - 1] + 1; r <= Fidx[newv[cnt - 1]]; r++ )
            {
              if ( dvadj[r - 1] == fc )
              {
                f1[f1cnt] = r - Fidx[newv[cnt - 1] - 1];
                f1cnt = f1cnt + 1;
              }
            }

            for ( r = 1; r <= f1cnt; r++ )
            {
              dvadj[Fidx[newv[cnt - 1] - 1] + f1[r - 1] - 1] = Nd;
            }
            dvadj[Fidx[fc - 1] + j - 1] = fc;
            free( f1 );
          }
        }

        curre = nexte;

        if ( Sixidx[curre - 1] == fin[k - 1] )
          nexte = Sixidx[curre - 1 + 2 * M_plan];
        else
          nexte = Sixidx[curre - 1 + 4 * M_plan];

        while ( bad_edges[nexte - 1] )
        {
          if ( Sixidx[nexte - 1] == fin[k - 1] )
            nexte = Sixidx[nexte - 1 + 2 * M_plan];
          else
            nexte = Sixidx[nexte - 1 + 4 * M_plan];
        }
      }

      Fidx[Nd] = Fidx[Nd - 1] + ( cnt );

      for ( i = 1; i <= cnt; i++ )
        dvadj[2 * Mdm + i - 1] = newv[i - 1];

      int totaldone = 1;
      bool* done = (bool*)malloc( Nd * sizeof( bool ) );
      int jj;
      for ( jj = 0; jj < Nd; jj++ )
        done[jj] = false;

      dist[k - 1 + ( Nd - 1 ) * d] = 0;
      done[Nd - 1] = true;
      int* looknext = (int*)malloc( ( Nd + 1 ) * sizeof( int ) );
      for ( jj = 0; jj < Nd + 1; jj++ )
        looknext[jj] = 0;
      looknext[0] = Nd;
      int looknextcount = 1;
      int index = 0;
      int count = 0;
      int target = current_crossings_best - requiredimprovement;
      while ( looknext[index] != 0 )
      {

        int v = looknext[index];
        if ( dist[k - 1 + ( v - 1 ) * d] >= target )
        {
          break;
        }

        for ( i = Fidx[v - 1]; i < Fidx[v]; i++ )
        {
          int check = ( dvadj[i] );
          if ( done[check - 1] == false )
          {
            looknext[looknextcount] = check;
            looknextcount = looknextcount + 1;
            dist[k - 1 + ( check - 1 ) * d] = dist[k - 1 + ( v - 1 ) * d] + 1;
            prev[k - 1 + ( check - 1 ) * d] = v;
            done[check - 1] = true;
            totaldone = totaldone + 1;
          }
        }
        index = index + 1;
      }

      free( dvadj );
      free( Fidx );
      free( newv );
      free( done );
      free( looknext );
    }

    for ( i = 1; i <= d; i++ )
    {
      for ( j = 1; j <= mcnt; j++ )
      {
        dist[i - 1 + ( merged_face_labels[j - 1] - 1 ) * d] = Mdm;
      }
      dist[i - 1 + ( Nd - 1 ) * d] = Mdm;
    }

    int* distsum = (int*)malloc( Nd * sizeof( int ) );
    for ( i = 1; i <= Nd; i++ )
      distsum[i - 1] = 0;

    for ( i = 1; i <= d; i++ )
      for ( j = 1; j <= Nd; j++ )
        distsum[j - 1] = distsum[j - 1] + dist[i - 1 + ( j - 1 ) * d];

    newcross = Mdm;
    for ( i = 1; i <= Nd; i++ )
      if ( distsum[i - 1] < newcross )
      {
        newcross = distsum[i - 1];
        newface = i;
      }
    if ( newcross + requiredimprovement < current_crossings_best )
    {
      improvementfound = 1;

      bool* checked = (bool*)malloc( Nd * sizeof( bool ) );
      for ( i = 0; i < Nd; i++ )
        checked[i] = false;
      maxpdist = 0;
      for ( i = 1; i <= d; i++ )
      {
        pdist[i - 1] = dist[i - 1 + ( newface - 1 ) * d];
        if ( pdist[i - 1] > maxpdist )
          maxpdist = pdist[i - 1];
      }

      pidx[0] = 0;
      pathsidx[0] = 0;
      for ( i = 1; i <= d; i++ )
      {
        pidx[i] = pidx[i - 1] + pdist[i - 1];
        pathsidx[i] = pathsidx[i - 1] + pdist[i - 1] + 1;
      }
      *paths = (int*)malloc( pathsidx[d] * sizeof( int ) );
      for ( i = 0; i < pathsidx[d]; i++ )
        ( *paths )[i] = 0;
      if ( maxpdist > 0 )
      {
        *pathedge = (int*)malloc( pidx[d] * sizeof( int ) );
        for ( i = 0; i < pidx[d]; i++ )
          ( *pathedge )[i] = 0;
      }
      else
      {
        *pathedge = (int*)malloc( d * sizeof( int ) );
        for ( i = 0; i < d; i++ )
          ( *pathedge )[i] = 0;
      }
      for ( k = 1; k <= d; k++ )
      {
        int nn = newface;
        ( *paths )[pathsidx[k - 1]] = nn;

        for ( j = 1; j <= pdist[k - 1] - 1; j++ )
        {
          nn = prev[k - 1 + ( nn - 1 ) * d];
          ( *paths )[pathsidx[k - 1] + j] = nn;
          if ( checked[nn - 1] )
            treecheck = 1;
          checked[nn - 1] = true;
        }

        if ( pdist[k - 1] > 0 )
        {
          int r;
          for ( r = 1; r <= maxdegfin; r++ )
          {
            if ( delv[k - 1 + ( r - 1 ) * d] == 0 )
              break;

            bool breakloop = false;

            for ( i = Fidxm[delv[k - 1 + ( r - 1 ) * d] - 1] + 1; i <= Fidxm[delv[k - 1 + ( r - 1 ) * d]]; i++ )
              if ( dvadjm[i - 1] == ( *paths )[pathsidx[k - 1] + ( pdist[k - 1] - 1 )] )
              {
                ( *paths )[pathsidx[k - 1] + ( pdist[k - 1] )] = delv[k - 1 + ( r - 1 ) * d];
                if ( checked[delv[k - 1 + ( r - 1 ) * d] - 1] )
                  treecheck = 1;
                checked[delv[k - 1 + ( r - 1 ) * d] - 1] = true;
                ( *pathedge )[pidx[k - 1] + ( pdist[k - 1] - 1 )] = deadjm[i - 1];
                pathedgeonce[deadjm[i - 1] - 1] = k;
                pathedgeonce[deadjm[i - 1] - 1 + maxdeadjm] = pdist[k - 1];
                breakloop = true;
                break;
              }
            if ( breakloop )
              break;
          }
        }
      }

      if ( maxpdist > 0 )
      {
        for ( k = 1; k <= d; k++ )
        {
          for ( j = 1; j <= pdist[k - 1] - 1; j++ )
          {
            for ( i = Fidxm[( *paths )[pathsidx[k - 1] + ( j - 1 )] - 1] + 1; i <= Fidxm[( *paths )[pathsidx[k - 1] + ( j - 1 )]]; i++ )
            {
              if ( dvadjm[i - 1] == ( *paths )[pathsidx[k - 1] + j] )
              {
                ( *pathedge )[pidx[k - 1] + ( j - 1 )] = deadjm[i - 1];
                pathedgeonce[deadjm[i - 1] - 1] = k;
                pathedgeonce[deadjm[i - 1] - 1 + maxdeadjm] = j;
                break;
              }
            }
          }
        }
      }

      free( checked );
    }
    else
    {
      *paths = (int*)malloc( d * sizeof( int ) );
      *pathedge = (int*)malloc( d * sizeof( int ) );
      for ( i = 1; i <= d; i++ )
      {
        ( *paths )[i - 1] = 0;
        ( *pathedge )[i - 1] = 0;
      }
    }
    free( dist );
    free( prev );
    free( distsum );
  }
  else
  {
    used_bigface = 1;
  }

  int maxmaxpaths = 0;
  for ( i = 1; i <= pathsidx[d]; i++ )
    if ( ( *paths )[i - 1] > maxmaxpaths )
      maxmaxpaths = ( *paths )[i - 1];
  outputs[0] = newcross;
  outputs[1] = treecheck;
  outputs[2] = newface;
  outputs[3] = maxpdist;
  outputs[4] = improvementfound;
  outputs[5] = used_bigface;
  outputs[6] = maxmaxpaths;
  outputs[7] = bigface_fails;
  outputs[8] = scheme;

  return bigface_success;
}

void Merge_Dual( int Nd, int* Fidx, int* dualel, int M, bool* bad_edge, int* dvadj, int* ddeg, int* deadj, int num_bad_edges, int* bad_edge_labels, int* newd, int* newfidx, int* newel, int* pooledgemap, int* merged_face_labels, int* newe, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Merge_Dual                                                             */
  /*                                                                         */
  /*  For each vertex considered for moving, we imagine that the vertex is   */
  /*  removed (along with its incident edges), and then a the planarised     */
  /*  graph is found and its faces identified. However, this is inefficient. */
  /*  Since we only really care about the dual graph (so we can compute      */
  /*  shortest paths), we instead just planarise the original graph, find    */
  /*  its faces, and hence the dual graph. Then for each vertex, we simply   */
  /*  update the dual graph. The removal of edges corresponds to the merging */
  /*  of faces, and hence the merging of vertices in the dual graph. Later,  */
  /*  we map back to the original faces by running Find_Merge_Faces.         */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;

  for ( i = 0; i < Fidx[Nd]; i++ )
  {
    newd[i] = 0;
    newe[i] = 0;
  }

  for ( i = 0; i < Nd + 1 + num_bad_edges; i++ )
    newfidx[0] = 0;

  for ( i = 0; i < 2 * M; i++ )
    newel[i] = 0;

  int* rel_chk = (int*)malloc( Nd * sizeof( int ) );
  for ( i = 0; i < Nd; i++ )
    rel_chk[i] = 0;

  int poolcnt = 0;

  for ( i = 0; i < 2 * num_bad_edges; i++ )
  {
    pooledgemap[i] = 0;
    merged_face_labels[i] = 0;
  }

  int mcnt = 0;

  int* degpool = (int*)malloc( Nd * sizeof( int ) );
  for ( i = 0; i < Nd; i++ )
    degpool[i] = ddeg[i];

  bool* face_checked = (bool*)malloc( Nd * sizeof( bool ) );
  for ( i = 0; i < Nd; i++ )
  {
    face_checked[i] = false;
  }

  for ( j = 1; j <= num_bad_edges; j++ )
  {
    i = bad_edge_labels[j - 1];
    degpool[dualel[i - 1] - 1] = degpool[dualel[i - 1] - 1] - 1;
    degpool[dualel[i - 1 + M] - 1] = degpool[dualel[i - 1 + M] - 1] - 1;
    int q1 = rel_chk[dualel[i - 1] - 1];
    int q2 = rel_chk[dualel[i - 1 + M] - 1];
    if ( q1 > 0 && q2 > 0 && q1 != q2 )
    {
      if ( q1 > q2 )
      {
        int qt = q1;
        q1 = q2;
        q2 = qt;
      }
      for ( k = 1; k <= Nd; k++ )
        if ( rel_chk[k - 1] == q2 )
          rel_chk[k - 1] = q1;
        else
        {
          if ( rel_chk[k - 1] > q2 )
            rel_chk[k - 1]--;
        }
      for ( k = q2; k <= poolcnt - 1; k++ )
        pooledgemap[k - 1] = pooledgemap[k];
      pooledgemap[poolcnt - 1] = 0;
      poolcnt--;
    }
    else
    {
      if ( q1 == 0 )
      {
        if ( q2 == 0 )
        {
          poolcnt++;
          rel_chk[dualel[i - 1] - 1] = poolcnt;
          rel_chk[dualel[i - 1 + M] - 1] = poolcnt;
          pooledgemap[poolcnt - 1] = dualel[i - 1];
          q2 = poolcnt;
          if ( face_checked[dualel[i - 1 + M] - 1] == false )
          {
            mcnt++;
            merged_face_labels[mcnt - 1] = dualel[i - 1 + M];
            face_checked[dualel[i - 1 + M] - 1] = true;
          }
        }
        else
        {
          rel_chk[dualel[i - 1] - 1] = q2;
        }
        q1 = poolcnt;
        if ( face_checked[dualel[i - 1] - 1] == false )
        {
          mcnt++;
          merged_face_labels[mcnt - 1] = dualel[i - 1];
          face_checked[dualel[i - 1] - 1] = true;
        }
      }
      if ( q2 == 0 )
      {
        rel_chk[dualel[i - 1 + M] - 1] = q1;
        if ( face_checked[dualel[i - 1 + M] - 1] == false )
        {
          mcnt++;
          merged_face_labels[mcnt - 1] = dualel[i - 1 + M];
          face_checked[dualel[i - 1 + M] - 1] = true;
        }
      }
    }
  }

  int* pooldeg = (int*)malloc( poolcnt * sizeof( int ) );
  for ( i = 0; i < poolcnt; i++ )
    pooldeg[i] = 0;

  int maxFsize = 0;
  int maxFindex = -1;
  int sumpooldeg = 0;

  for ( i = 1; i <= mcnt; i++ )
  {
    pooldeg[rel_chk[merged_face_labels[i - 1] - 1] - 1] = pooldeg[rel_chk[merged_face_labels[i - 1] - 1] - 1] + degpool[merged_face_labels[i - 1] - 1];
    sumpooldeg = sumpooldeg + degpool[merged_face_labels[i - 1] - 1];
    if ( pooldeg[rel_chk[merged_face_labels[i - 1] - 1] - 1] == maxFsize )
    {
      if ( maxFindex > rel_chk[merged_face_labels[i - 1] - 1] )
        maxFindex = rel_chk[merged_face_labels[i - 1] - 1];
    }
    if ( pooldeg[rel_chk[merged_face_labels[i - 1] - 1] - 1] > maxFsize )
    {
      maxFsize = pooldeg[rel_chk[merged_face_labels[i - 1] - 1] - 1];
      maxFindex = rel_chk[merged_face_labels[i - 1] - 1];
    }
  }

  memcpy( newel, dualel, 2 * M * sizeof( int ) );

  int* poole = (int*)malloc( sumpooldeg * sizeof( int ) );
  int* poolf = (int*)malloc( sumpooldeg * sizeof( int ) );
  for ( i = 0; i < sumpooldeg; i++ )
  {
    poole[i] = 0;
    poolf[i] = 0;
  }

  int* poolidx = (int*)malloc( ( poolcnt + 1 ) * sizeof( int ) );
  poolidx[0] = 0;
  for ( i = 2; i <= poolcnt + 1; i++ )
    poolidx[i - 1] = poolidx[i - 2] + pooldeg[i - 2];

  int* poolcurr = (int*)malloc( poolcnt * sizeof( int ) );
  for ( i = 0; i < poolcnt; i++ )
    poolcurr[i] = 0;

  for ( i = 1; i <= M; i++ )
  {
    if ( rel_chk[dualel[i - 1] - 1] > 0 && bad_edge[i - 1] == false )
    {
      int q1 = rel_chk[dualel[i - 1] - 1];
      newel[i - 1] = Nd + q1;
      poolcurr[q1 - 1]++;
      if ( rel_chk[dualel[i - 1 + M] - 1] == 0 )
        poolf[poolidx[q1 - 1] + poolcurr[q1 - 1] - 1] = dualel[i - 1 + M];
      else
        poolf[poolidx[q1 - 1] + poolcurr[q1 - 1] - 1] = Nd + rel_chk[dualel[i - 1 + M] - 1];

      poole[poolidx[q1 - 1] + poolcurr[q1 - 1] - 1] = i;
    }
    if ( rel_chk[dualel[i - 1 + M] - 1] > 0 && bad_edge[i - 1] == false )
    {
      int q2 = rel_chk[dualel[i - 1 + M] - 1];
      newel[i - 1 + M] = Nd + q2;
      poolcurr[q2 - 1]++;
      if ( rel_chk[dualel[i - 1] - 1] == 0 )
        poolf[poolidx[q2 - 1] + poolcurr[q2 - 1] - 1] = dualel[i - 1];
      else
        poolf[poolidx[q2 - 1] + poolcurr[q2 - 1] - 1] = Nd + rel_chk[dualel[i - 1] - 1];
      poole[poolidx[q2 - 1] + poolcurr[q2 - 1] - 1] = i;
    }
  }

  int cnt = 0;
  int cntstart = 0;
  int istart = 0;
  int newmaxFindex = -1;
  for ( i = 1; i <= Nd; i++ )
  {
    if ( rel_chk[i - 1] == 0 )
    {
      if ( istart == 0 )
      {
        cntstart = cnt;
        istart = i;
      }
      cnt = cnt + Fidx[i] - Fidx[i - 1];
      if ( Fidx[i] - Fidx[i - 1] > maxFsize )
      {
        maxFsize = Fidx[i] - Fidx[i - 1];
        newmaxFindex = i;
      }
    }
    else
    {
      if ( istart > 0 )
      {
        for ( j = 1; j <= cnt - cntstart; j++ )
        {
          newd[cntstart + j - 1] = dvadj[Fidx[istart - 1] + j - 1];
          newe[cntstart + j - 1] = deadj[Fidx[istart - 1] + j - 1];
        }
        istart = 0;
      }
    }
    newfidx[i] = cnt;
  }

  if ( newmaxFindex < 0 )
    maxFindex = maxFindex + Nd;
  else
    maxFindex = newmaxFindex;

  if ( rel_chk[Nd - 1] == 0 )
  {
    for ( i = 1; i <= cnt - cntstart; i++ )
    {
      newd[cntstart + i - 1] = dvadj[Fidx[istart - 1] + i - 1];
      newe[cntstart + i - 1] = deadj[Fidx[istart - 1] + i - 1];
    }
  }

  for ( i = 1; i <= cnt; i++ )
  {
    if ( rel_chk[newd[i - 1] - 1] > 0 )
      newd[i - 1] = Nd + rel_chk[newd[i - 1] - 1];
  }

  for ( i = 1; i <= sumpooldeg; i++ )
  {
    newe[cnt + i - 1] = poole[i - 1];
    newd[cnt + i - 1] = poolf[i - 1];
  }

  for ( i = 1; i <= poolcnt; i++ )
    newfidx[Nd + i] = poolidx[i] + cnt;

  int maxnewe = 0;
  for ( i = 1; i <= cnt + sumpooldeg; i++ )
    if ( newe[i] > maxnewe )
      maxnewe = newe[i];

  outputs[0] = Nd + poolcnt;
  outputs[1] = newfidx[Nd + poolcnt] / 2;
  outputs[2] = maxFindex;
  outputs[3] = poolcnt;
  outputs[4] = mcnt;
  outputs[5] = cnt;
  outputs[6] = sumpooldeg;
  outputs[7] = maxnewe;

  free( degpool );
  free( face_checked );
  free( pooldeg );
  free( poole );
  free( poolf );
  free( poolidx );
  free( poolcurr );
  free( rel_chk );
}

void Prep_Vars( int* deg, int* incidente, int* Sixidx, int* elist, int* eidx, int* Clabel, int* Cidx, int M, int M_plan, int vm, int current_crossings, int** fin, int** fi, int** bad_edge_labels, bool** bad_edges, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Prep_Vars                                                              */
  /*                                                                         */
  /*  A number of variables need to be prepared so that Merge_Dual can be    */
  /*  executed, which is done in this subroutine.                            */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;

  int d = deg[vm - 1];
  *fin = (int*)malloc( d * sizeof( int ) );
  *fi = (int*)malloc( d * sizeof( int ) );
  for ( j = 1; j <= d; j++ )
  {
    ( *fin )[j - 1] = 0;
    ( *fi )[j - 1] = 0;
  }

  ( *fi )[0] = incidente[vm - 1];

  if ( Sixidx[( *fi )[0] - 1] == vm )
    ( *fin )[0] = Sixidx[( *fi )[0] - 1 + M];
  else
    ( *fin )[0] = Sixidx[( *fi )[0] - 1];

  for ( j = 2; j <= d; j++ )
  {
    if ( Sixidx[( *fi )[j - 2] - 1] == vm )
      ( *fi )[j - 1] = Sixidx[( *fi )[j - 2] - 1 + 2 * M];
    else
      ( *fi )[j - 1] = Sixidx[( *fi )[j - 2] - 1 + 4 * M];

    if ( Sixidx[( *fi )[j - 1] - 1] == vm )
      ( *fin )[j - 1] = Sixidx[( *fi )[j - 1] - 1 + M];
    else
      ( *fin )[j - 1] = Sixidx[( *fi )[j - 1] - 1];
  }

  int maxdegfin = 0;
  for ( j = 1; j <= d; j++ )
    if ( deg[( *fin )[j - 1] - 1] > maxdegfin )
      maxdegfin = deg[( *fin )[j - 1] - 1];

  int bad_edge_labels_size = 0;
  for ( j = 1; j <= d; j++ )
    bad_edge_labels_size = bad_edge_labels_size + eidx[( *fi )[j - 1]] - eidx[( *fi )[j - 1] - 1];

  *bad_edge_labels = (int*)malloc( bad_edge_labels_size * sizeof( int ) );
  for ( j = 1; j <= bad_edge_labels_size; j++ )
    ( *bad_edge_labels )[j - 1] = 0;

  *bad_edges = (bool*)malloc( M_plan * sizeof( bool ) );
  for ( j = 1; j <= M_plan; j++ )
    ( *bad_edges )[j - 1] = false;

  int cnt1 = 0;

  for ( j = 1; j <= d; j++ )
  {
    for ( k = eidx[( *fi )[j - 1] - 1] + 1; k <= eidx[( *fi )[j - 1]]; k++ )
    {
      ( *bad_edges )[elist[k - 1] - 1] = true;
      ( *bad_edge_labels )[cnt1 + k - eidx[( *fi )[j - 1] - 1] - 1] = elist[k - 1];
    }
    cnt1 = cnt1 + eidx[( *fi )[j - 1]] - eidx[( *fi )[j - 1] - 1];
  }

  int requiredimprovement = current_crossings;
  for ( i = 1; i <= d; i++ )
    for ( j = Cidx[( *fi )[i - 1] - 1] + 1; j <= Cidx[( *fi )[i - 1]]; j++ )
    {
      bool foundClabelj = false;
      for ( k = 1; k <= d; k++ )
        if ( ( *fi )[k - 1] == Clabel[j - 1] )
        {
          foundClabelj = true;
          if ( ( *fi )[i - 1] > Clabel[j - 1] )
            requiredimprovement--;
          break;
        }
      if ( foundClabelj == false )
        requiredimprovement--;
    }

  outputs[0] = d;
  outputs[1] = requiredimprovement;
  outputs[2] = bad_edge_labels_size;
  outputs[3] = maxdegfin;
}

int Faces( int* Sixidx, int M, int N, int* fi, int fi_size, int* F, int* Fidx, int* Fdir, int* dvadj, int* ddeg, int* dualel, int* deadj )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Faces                                                                  */
  /*                                                                         */
  /*  Once the graph has been planarised, we then find its faces, and hence  */
  /*  the dual graph. Although we technically need this for graphs with a    */
  /*  vertex (and its incident edges) removed, in practice we only do it for */
  /*  the original graph at the start of each iteration, and then perform    */
  /*  updates on the dual graph for each removed vertex using Merge_Dual.    */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  for ( i = 0; i < 2 * M; i++ )
  {
    dualel[i] = 0;
    F[i] = 0;
    Fdir[i] = 0;
    dvadj[i] = 0;
    deadj[i] = 0;
  }

  for ( i = 0; i < ( 2 * M + 5 ) / 3; i++ )
    Fidx[0] = 0;

  for ( i = 0; i < 2 + M - N; i++ )
    ddeg[i] = 0;

  int fnum = 0;
  int fcount = 0;

  int* queue = (int*)malloc( 2 * M * sizeof( int ) );
  queue[0] = 1;
  for ( i = 1; i < 2 * M; i++ )
    queue[i] = 0;

  bool increase_queue1 = false;
  for ( i = 0; i < fi_size; i++ )
    if ( fi[i] == queue[0] )
    {
      increase_queue1 = true;
      break;
    }

  while ( increase_queue1 )
  {
    queue[0]++;
    increase_queue1 = false;
    for ( i = 0; i < fi_size; i++ )
      if ( fi[i] == queue[0] )
      {
        increase_queue1 = true;
        break;
      }
  }

  int* considered = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
    considered[i] = 0;

  int qi = 0;
  int qni = 1;
  while ( true )
  {
    int startv = -1;
    int starte = -1;
    int currd = 0;
    int currv = -1;
    int curre = -1;

    qi++;
    if ( qi > qni )
      break;

    int ed = queue[qi - 1];
    if ( ed > 0 )
    {
      starte = ed;
      if ( considered[starte - 1] == 2 )
        continue;

      fnum++;
      startv = Sixidx[ed - 1];
      dualel[ed - 1] = fnum;
      currd = 1;
      if ( considered[starte - 1] == 0 )
      {
        qni++;
        queue[qni - 1] = -ed;
      }
      considered[starte - 1] = 1 - considered[starte - 1];
    }
    else
    {
      starte = -ed;
      if ( considered[starte - 1] == 2 )
        continue;

      fnum++;
      startv = Sixidx[-ed - 1 + M];
      dualel[-ed - 1 + M] = fnum;
      currd = -1;
      if ( considered[starte - 1] == 0 )
      {
        qni++;
        queue[qni - 1] = -ed;
      }
      considered[starte - 1] = -1 + 3 * considered[starte - 1];
    }
    curre = starte;
    currv = startv;
    fcount++;
    F[fcount - 1] = curre;
    Fdir[fcount - 1] = currd;
    while ( true )
    {
      if ( currd == 1 )
      {
        currv = Sixidx[curre - 1 + M];
        curre = Sixidx[curre - 1 + 5 * M];
      }
      else
      {
        currv = Sixidx[curre - 1];
        curre = Sixidx[curre - 1 + 3 * M];
      }

      if ( Sixidx[curre - 1] == currv )
        currd = 1;
      else
        currd = -1;

      if ( currv == startv && curre == starte )
      {
        Fidx[fnum] = fcount;
        break;
      }
      else
      {
        if ( considered[curre - 1] == -currd )
          considered[curre - 1] = 2;
        else
        {
          considered[curre - 1] = currd;
          qni++;
          queue[qni - 1] = -currd * curre;
        }

        if ( currd == 1 )
          dualel[curre - 1] = fnum;
        else
          dualel[curre - 1 + M] = fnum;

        fcount++;

        F[fcount - 1] = curre;
        Fdir[fcount - 1] = currd;
      }
    }
  }

  for ( i = 1; i <= M; i++ )
  {
    bool skip = false;
    for ( j = 1; j <= fi_size; j++ )
      if ( fi[j - 1] == i )
      {
        skip = true;
        break;
      }
    if ( skip )
      continue;

    ddeg[dualel[i - 1] - 1]++;
    dvadj[Fidx[dualel[i - 1] - 1] + ddeg[dualel[i - 1] - 1] - 1] = dualel[i - 1 + M];
    deadj[Fidx[dualel[i - 1] - 1] + ddeg[dualel[i - 1] - 1] - 1] = i;

    ddeg[dualel[i - 1 + M] - 1]++;
    dvadj[Fidx[dualel[i - 1 + M] - 1] + ddeg[dualel[i - 1 + M] - 1] - 1] = dualel[i - 1];
    deadj[Fidx[dualel[i - 1 + M] - 1] + ddeg[dualel[i - 1 + M] - 1] - 1] = i;
  }

  free( considered );
  free( queue );
  return fnum;
}

void Planarise_Sixidx( int* Sixidx, int* Clabel, int* Cidx, int* Dlabel, int N, int M, int* Sixidx_plan, int* planvlabels, int* elist, int* eidx, int* origedgelabel, int* origedgeorder, int* incidente_plan, int* newdim )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Planarise_Sixidx                                                       */
  /*                                                                         */
  /*  The graph is planarised, and the combinatorial embedding for this new  */
  /*  planar graph is produced.                                              */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int k;

  int* otherplanvlabel = (int*)malloc( Cidx[M] * sizeof( int ) );

  for ( i = 0; i < Cidx[M]; i++ )
  {
    planvlabels[i] = 0;
    otherplanvlabel[i] = 0;
  }
  for ( i = 0; i < M + Cidx[M]; i++ )
    elist[i] = 0;

  for ( i = 0; i < M + 1; i++ )
    eidx[i] = i + Cidx[i];

  for ( i = 0; i < M; i++ )
    elist[eidx[i]] = i + 1;

  int oldN = N;
  int oldM = M;

  for ( i = 0; i < M; i++ )
  {
    origedgelabel[i] = i + 1;
    origedgeorder[i] = 1;
  }
  for ( i = 0; i < Cidx[M]; i++ )
  {
    origedgelabel[M + i] = 0;
    origedgeorder[M + i] = 0;
  }

  for ( i = 0; i < N; i++ )
    incidente_plan[i] = 0;

  int nss = M + Cidx[M];

  for ( i = 0; i < M; i++ )
    for ( j = 0; j < 6; j++ )
      Sixidx_plan[i + j * nss] = Sixidx[i + j * M];

  for ( i = 0; i < Cidx[M]; i++ )
    for ( j = 0; j < 6; j++ )
      Sixidx_plan[M + i + j * nss] = 0;

  for ( i = 1; i <= oldM; i++ )
  {
    int tt = -1;
    if ( Cidx[i] - Cidx[i - 1] != 0 )
      tt = Sixidx_plan[( i - 1 ) + nss];

    for ( j = 1; j <= Cidx[i] - Cidx[i - 1]; j++ )
    {
      M = M + 1;
      elist[eidx[i - 1] + j] = M;
      origedgelabel[M - 1] = i;
      origedgeorder[M - 1] = j + 1;

      if ( planvlabels[Cidx[i - 1] + j - 1] == 0 )
      {
        N = N + 1;
        planvlabels[Cidx[i - 1] + j - 1] = N;
        int ce = Clabel[Cidx[i - 1] + j - 1];

        for ( k = Cidx[ce - 1] + 1; k <= Cidx[ce]; k++ )
        {
          if ( Clabel[k - 1] == i )
          {
            planvlabels[k - 1] = N;
            otherplanvlabel[Cidx[i - 1] + j - 1] = k;
            otherplanvlabel[k - 1] = Cidx[i - 1] + j;
          }
        }
      }
    }

    if ( Cidx[i] - Cidx[i - 1] > 0 )
    {
      Sixidx_plan[M - 1] = tt;
      Sixidx_plan[M - 1 + nss] = planvlabels[Cidx[i] - 1];
      Sixidx_plan[M - 1 + 2 * nss] = Sixidx_plan[i - 1 + 4 * nss];
      Sixidx_plan[M - 1 + 3 * nss] = Sixidx_plan[i - 1 + 5 * nss];
      Sixidx_plan[M - 1 + 4 * nss] = 0;
      Sixidx_plan[M - 1 + 5 * nss] = 0;
      Sixidx_plan[i - 1 + nss] = planvlabels[Cidx[i - 1]];
    }

    for ( j = 1; j <= Cidx[i] - Cidx[i - 1] - 1; j++ )
    {
      if ( planvlabels[Cidx[i - 1] + j - 1] < planvlabels[Cidx[i - 1] + j] )
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1] = planvlabels[Cidx[i - 1] + j - 1];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + nss] = planvlabels[Cidx[i - 1] + j];
      }
      else
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1] = planvlabels[Cidx[i - 1] + j];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + nss] = planvlabels[Cidx[i - 1] + j - 1];
      }
    }
  }

  for ( i = 1; i <= M; i++ )
  {
    if ( Sixidx_plan[i - 1] <= oldN )
    {
      int ec = Sixidx_plan[i - 1 + 2 * nss];
      if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1] )
        Sixidx_plan[i - 1 + 2 * nss] = elist[eidx[ec] - 1];

      ec = Sixidx_plan[i - 1 + 3 * nss];
      if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1] )
        Sixidx_plan[i - 1 + 3 * nss] = elist[eidx[ec] - 1];
    }
    if ( Sixidx_plan[i - 1 + nss] <= oldN )
    {
      int ec = Sixidx_plan[i - 1 + 4 * nss];
      if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1 + nss] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1 + nss] )
        Sixidx_plan[i - 1 + 4 * nss] = elist[eidx[ec] - 1];

      ec = Sixidx_plan[i - 1 + 5 * nss];
      if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1 + nss] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1 + nss] )
        Sixidx_plan[i - 1 + 5 * nss] = elist[eidx[ec] - 1];
    }
  }

  for ( i = 1; i <= oldM; i++ )
  {
    for ( j = 1; j <= Cidx[i] - Cidx[i - 1]; j++ )
    {
      int othere = Clabel[Cidx[i - 1] + j - 1];
      int f1 = otherplanvlabel[Cidx[i - 1] + j - 1] - Cidx[othere - 1];
      int a1 = 1;
      int a2 = 0;

      if ( Dlabel[Cidx[i - 1] + j - 1] != 1 )
      {
        a1 = 0;
        a2 = 1;
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1] )
      {
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 2 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 3 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
      }
      else
      {
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 4 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 5 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[i - 1] + j] - 1] )
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 2 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 3 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
      }
      else
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 4 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 5 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1] )
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 2 * nss] = elist[eidx[i - 1] + j - 1 + a2];
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 3 * nss] = elist[eidx[i - 1] + j - 1 + a1];
      }
      else
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 4 * nss] = elist[eidx[i - 1] + j - 1 + a2];
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 5 * nss] = elist[eidx[i - 1] + j - 1 + a1];
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[othere - 1] + f1] - 1] )
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 2 * nss] = elist[eidx[i - 1] + j - 1 + a1];
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 3 * nss] = elist[eidx[i - 1] + j - 1 + a2];
      }
      else
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 4 * nss] = elist[eidx[i - 1] + j - 1 + a1];
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 5 * nss] = elist[eidx[i - 1] + j - 1 + a2];
      }
    }
  }

  for ( i = 1; i <= M; i++ )
  {
    if ( Sixidx_plan[i - 1] <= oldN && incidente_plan[Sixidx_plan[i - 1] - 1] == 0 )
      incidente_plan[Sixidx_plan[i - 1] - 1] = i;
    if ( Sixidx_plan[i - 1 + nss] <= oldN && incidente_plan[Sixidx_plan[i - 1 + nss] - 1] == 0 )
      incidente_plan[Sixidx_plan[i - 1 + nss] - 1] = i;
  }

  newdim[0] = N;
  newdim[1] = M;

  free( otherplanvlabel );
}

void Quick_Cross_Main_Loop( int** incidente, int** Sixidx, int** Clabel, int** Dlabel, int** Cidx, int N, int M, int** deg, int current_crossings, int min_scheme, int bigface_depth, int stop_crossings, bool ignore_stop_crossings, int orig_N, int* outputs2, int verbose )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Quick_Cross_Main_Loop                                                  */
  /*                                                                         */
  /*  This is the main loop which performs iterations of QuickCross. It      */
  /*  functions as follows:                                                  */
  /*                                                                         */
  /*  Planarise the graph in its current embedding (Planarise_Sixidx)        */
  /*  Find the faces, and hence dual graph of the planarised graph (Faces)   */
  /*  Then, consider each vertex for moving, and:                            */
  /*    Find the merged dual graph (Prep_Vars and Merge_Dual)                */
  /*    Find the best (or an improved) destination face (Get_Distances)      */
  /*    If min_scheme = 1 or 2, break loop if an improvement is found.       */
  /*    If min_scheme = 3 finish loop and determine the biggest improvement. */
  /*  If no improvement was found, STOP. Otherwise:                          */
  /*    Relocate vertex and update graph data (Find_Merge_Face, New_Clabel)  */
  /*    Check incident edge crosses or if subdivisions are needed (Subd_Try) */
  /*    Check to see if previous subdivisions can be removed (Undo_Sub)      */
  /*                                                                         */
  /*  Each iteration should take at most O(kM) time where M is the number of */
  /*  edges in the original graph, and k is the current number of crossings. */
  /*                                                                         */
  /***************************************************************************/

  if ( verbose > 0 )
    fprintf( stdout, "\nStarting Main Loop with %d crossings:\n\n", current_crossings );
  int debug = 0;
  int i;
  int j;
  int k;

  int iters = 0;
  int vm_best = 0;
  int bigface_fails = 0;

  bool breakloop = false;

  while ( true )
  {
    if ( ( ignore_stop_crossings == false && current_crossings <= stop_crossings ) || current_crossings == 0 )
    {
      if ( verbose > 0 )
        fprintf( stdout, "Final crossings for this component: %d\n\n", current_crossings );
      break;
    }

    int current_crossings_best = current_crossings;
    iters++;
    int improvementfound = 0;

    if ( verbose > 1 )
      fprintf( stdout, "\nStarting iteration %d:", iters );

    if ( debug == 1 )
      fprintf( stdout, "Running Planarise_Sixidx\n" );
    int* Sixidx_plan = (int*)malloc( ( 6 * M + 6 * ( *Cidx )[M] ) * sizeof( int ) );
    int* planvlabels = (int*)malloc( ( ( *Cidx )[M] ) * sizeof( int ) );
    int* elist = (int*)malloc( ( M + ( *Cidx )[M] ) * sizeof( int ) );
    int* eidx = (int*)malloc( ( M + 1 ) * sizeof( int ) );
    int* origedgelabel = (int*)malloc( ( M + ( *Cidx )[M] ) * sizeof( int ) );
    int* origedgeorder = (int*)malloc( ( M + ( *Cidx )[M] ) * sizeof( int ) );
    int* incidente_plan = (int*)malloc( ( N ) * sizeof( int ) );
    int newdim[2];
    Planarise_Sixidx( *Sixidx, *Clabel, *Cidx, *Dlabel, N, M, Sixidx_plan, planvlabels, elist, eidx, origedgelabel, origedgeorder, incidente_plan, newdim );
    if ( debug == 1 )
      fprintf( stdout, "Finished Running Planarise_Sixidx\n" );
    int N_plan = newdim[0];
    int M_plan = newdim[1];

    if ( debug == 1 )
      fprintf( stdout, "Running Faces\n" );
    int* F = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* Fidx = (int*)malloc( ( ( 2 * M_plan + 5 ) / 3 ) * sizeof( int ) );
    int* Fdir = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* dvadj = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* ddeg = (int*)malloc( ( 2 + M_plan - N_plan ) * sizeof( int ) );
    int* dualel = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* deadj = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int Nd = Faces( Sixidx_plan, M_plan, N_plan, NULL, 0, F, Fidx, Fdir, dvadj, ddeg, dualel, deadj );

    int* delv_best;
    int treecheck_best;
    int* pathedge_best;
    int* pidx_best;
    int* pathedgeonce_best;
    int* paths_best;
    int* pathsidx_best;
    int* pdist_best;
    int* fin_best;
    int* fi_best;
    bool* bad_edges_best;
    int* dualel_best;
    int newface_best;
    int d_best;
    int* pooledgemap_best;
    int maxpdist_best;
    int maxmaxpaths_best;
    int maxdeadjm_best;

    bool assigned_best = false;
    bool breakinnerloop = false;

    int vm = vm_best;

    int requiredimprovement = -1;

    for ( i = 1; i <= N; i++ )
    {
      vm++;
      if ( vm > N )
        vm = vm - N;

      if ( verbose > 1 )
        fprintf( stdout, "\n Trying vm = %d... ", vm );
      if ( debug == 1 )
        fprintf( stdout, "Running Prep_Vars... " );
      int* fin;
      int* fi;
      int* bad_edge_labels;
      bool* bad_edges;
      int prepoutputs[4];
      Prep_Vars( *deg, *incidente, *Sixidx, elist, eidx, *Clabel, *Cidx, M, M_plan, vm, current_crossings, &fin, &fi, &bad_edge_labels, &bad_edges, prepoutputs );
      int d = prepoutputs[0];
      int requiredimprovement = prepoutputs[1];
      int num_bad_edges = prepoutputs[2];
      int maxdegfin = prepoutputs[3];

      if ( debug == 1 )
        fprintf( stdout, "Running Merge_Dual... " );
      int* newdvadj = (int*)malloc( ( Fidx[Nd] ) * sizeof( int ) );
      int* newfidx = (int*)malloc( ( Nd + 1 + num_bad_edges ) * sizeof( int ) );
      int* newdualel = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
      int* pooledgemap = (int*)malloc( ( 2 * num_bad_edges ) * sizeof( int ) );
      int* merged_face_labels = (int*)malloc( ( 2 * num_bad_edges ) * sizeof( int ) );
      int* neweadj = (int*)malloc( ( Fidx[Nd] ) * sizeof( int ) );
      int mdoutputs[8];
      Merge_Dual( Nd, Fidx, dualel, M_plan, bad_edges, dvadj, ddeg, deadj, num_bad_edges, bad_edge_labels, newdvadj, newfidx, newdualel, pooledgemap, merged_face_labels, neweadj, mdoutputs );
      int newNd = mdoutputs[0];
      int newMd = mdoutputs[1];
      int maxFindex = mdoutputs[2];
      int poolcnt = mdoutputs[3];
      int mcnt = mdoutputs[4];
      int cnt = mdoutputs[5];
      int sumpooldeg = mdoutputs[6];
      int maxdeadjm = mdoutputs[7];

      if ( debug == 1 )
        fprintf( stdout, "Running Get_Distances...\n" );
      int* pathedge;
      int* paths;
      int* pathedgeonce = (int*)malloc( ( 2 * maxdeadjm ) * sizeof( int ) );
      int* pdist = (int*)malloc( d * sizeof( int ) );
      int* pidx = (int*)malloc( ( d + 1 ) * sizeof( int ) );
      int* pathsidx = (int*)malloc( ( d + 1 ) * sizeof( int ) );
      int* delv = (int*)malloc( d * maxdegfin * sizeof( int ) );
      int gdoutputs[9];
      int bigface_success = Get_Distances( newdvadj, newfidx, newNd, newMd, M_plan, Sixidx_plan, newdualel, requiredimprovement, d, fin, incidente_plan, *deg, bad_edges, current_crossings_best, neweadj, mcnt, merged_face_labels, maxFindex, min_scheme, bigface_depth, maxdegfin, maxdeadjm, &pathedge, pidx, pathedgeonce, &paths, pathsidx, pdist, delv, gdoutputs, bigface_fails );
      int newcross = gdoutputs[0];
      int treecheck = gdoutputs[1];
      int newface = gdoutputs[2];
      int maxpdist = gdoutputs[3];
      int tempimprovementfound = gdoutputs[4];
      int used_bigface = gdoutputs[5];
      int maxmaxpaths = gdoutputs[6];
      bigface_fails = gdoutputs[7];
      min_scheme = gdoutputs[8];

      if ( debug == 1 )
        fprintf( stdout, "Finished Get_Distances...\n" );

      if ( current_crossings - requiredimprovement - newcross > 0 )
      {
        if ( verbose > 1 )
        {
          if ( bigface_success >= 0 )
            fprintf( stdout, "Found improvement in bigface depth %d (%d crossings removed)", bigface_success, current_crossings - requiredimprovement - newcross );
          else
            fprintf( stdout, "Found improvement (%d crossings removed)", current_crossings - requiredimprovement - newcross );
        }
      }
      else
      {
        if ( verbose > 1 )
          fprintf( stdout, "No improvement possible." );
      }

      if ( requiredimprovement + newcross < current_crossings_best )
      {
        if ( verbose > 1 && min_scheme == 3 )
          fprintf( stdout, " (best so far)" );
        improvementfound = 1;
        if ( assigned_best )
        {
          free( delv_best );
          free( pathedge_best );
          free( pidx_best );
          free( pathedgeonce_best );
          free( paths_best );
          free( pathsidx_best );
          free( pdist_best );
          free( fin_best );
          free( fi_best );
          free( bad_edges_best );
          free( dualel_best );
          free( pooledgemap_best );
        }
        assigned_best = true;

        current_crossings_best = requiredimprovement + newcross;
        vm_best = vm;
        treecheck_best = treecheck;
        newface_best = newface;
        d_best = d;
        maxpdist_best = maxpdist;
        maxmaxpaths_best = maxmaxpaths;
        maxdeadjm_best = maxdeadjm;
        copy_array( delv, &delv_best, d * maxdegfin );
        if ( maxpdist == 0 )
          copy_array( pathedge, &pathedge_best, d );
        else
          copy_array( pathedge, &pathedge_best, pidx[d] );
        copy_array( pidx, &pidx_best, d + 1 );
        copy_array( pathedgeonce, &pathedgeonce_best, 2 * maxdeadjm );
        copy_array( paths, &paths_best, pathsidx[d] );
        copy_array( pathsidx, &pathsidx_best, d + 1 );
        copy_array( pdist, &pdist_best, d );
        copy_array( fin, &fin_best, d );
        copy_array( fi, &fi_best, d );
        copy_array_bool( bad_edges, &bad_edges_best, M_plan );
        copy_array( newdualel, &dualel_best, 2 * M_plan );
        copy_array( pooledgemap, &pooledgemap_best, 2 * num_bad_edges );

        if ( min_scheme < 3 )
          breakinnerloop = true;
      }

      if ( debug == 1 )
        fprintf( stdout, "Freeing Variables\n" );
      free( pathedge );
      free( paths );
      free( fin );
      free( fi );
      free( bad_edge_labels );
      free( bad_edges );
      free( newdvadj );
      free( newfidx );
      free( newdualel );
      free( pooledgemap );
      free( merged_face_labels );
      free( neweadj );
      free( pathedgeonce );
      free( pdist );
      free( pidx );
      free( pathsidx );
      free( delv );

      if ( breakinnerloop )
        break;
    }

    if ( improvementfound == 1 )
    {

      if ( debug == 1 )
        fprintf( stdout, "Running Find_Merge_Faces\n" );
      int* newf_best = (int*)malloc( Fidx[Nd] * sizeof( int ) );
      int* newfidx_best = (int*)malloc( Nd * sizeof( int ) );
      int* newdir_best = (int*)malloc( Fidx[Nd] * sizeof( int ) );
      int* rf_map_best = (int*)malloc( maxmaxpaths_best * sizeof( int ) );
      Find_Merge_Faces( Nd, M_plan, F, Fdir, Fidx, dualel_best, bad_edges_best, paths_best, pathsidx_best, pdist_best, d_best, Nd, pooledgemap_best, newf_best, newfidx_best, newdir_best, rf_map_best, maxmaxpaths_best );

      current_crossings = current_crossings_best;
      for ( i = 1; i <= d_best; i++ )
        for ( j = 1; j <= pdist_best[i - 1] + 1; j++ )
        {
          paths_best[pathsidx_best[i - 1] + ( j - 1 )] = rf_map_best[paths_best[pathsidx_best[i - 1] + ( j - 1 )] - 1];
        }

      newface_best = rf_map_best[newface_best - 1];
      if ( debug == 1 )
        fprintf( stdout, "Up to New Clabel\n" );

      New_Clabel( M, M_plan, Sixidx_plan, *Sixidx, *Clabel, *Dlabel, *Cidx, elist, eidx, pathedge_best, paths_best, pidx_best, pathsidx_best, pdist_best, maxpdist_best, d_best, fin_best, fi_best, newf_best, newdir_best, newfidx_best, vm_best, origedgelabel, origedgeorder, newface_best, pathedgeonce_best, treecheck_best, maxdeadjm_best );

      if ( verbose > 1 )
        fprintf( stdout, "\n" );
      if ( verbose > 0 )
        fprintf( stdout, "Iteration %d: Move vertex %d, crossings = %d\n", iters, vm_best, current_crossings );

      int stoutputs[3];

      Subd_Try( current_crossings, Sixidx, *Clabel, *Dlabel, Cidx, pathedge_best, N, M, origedgelabel, incidente, deg, pdist_best, fi_best, fin_best, vm_best, d_best, maxpdist_best, pidx_best, stoutputs, orig_N, false, verbose );
      N = stoutputs[0];
      M = stoutputs[1];
      current_crossings = stoutputs[2];
      free( newf_best );
      free( newfidx_best );
      free( newdir_best );
      free( rf_map_best );
    }
    else
    {
      if ( verbose > 1 )
        fprintf( stdout, "\n\n" );
      if ( verbose > 0 )
        fprintf( stdout, "Final crossings for this component: %d\n\n", current_crossings );
      breakloop = true;
    }

    free( Sixidx_plan );
    free( planvlabels );
    free( elist );
    free( eidx );
    free( origedgelabel );
    free( origedgeorder );
    free( incidente_plan );
    free( F );
    free( Fidx );
    free( Fdir );
    free( dvadj );
    free( ddeg );
    free( dualel );
    free( deadj );

    if ( assigned_best )
    {
      free( delv_best );
      free( pathedge_best );
      free( pidx_best );
      free( pathedgeonce_best );
      free( paths_best );
      free( pathsidx_best );
      free( pdist_best );
      free( fin_best );
      free( fi_best );
      free( bad_edges_best );
      free( dualel_best );
      free( pooledgemap_best );
    }

    if ( breakloop )
      break;
  }

  outputs2[0] = current_crossings;
  outputs2[1] = M;
  outputs2[2] = N;
}

void Subd_Try_PE( int current_crossings, int** Sixidx, int* Clabel, int* Dlabel, int** Cidx, int* pathedge, int N, int M, int* origedgelabel, int** incidente, int** deg, int* pdist, int* fi, int* fin, int vm, int d, int maxpdist, int* pidx, int* outputs, int** not_added_edges, int** checked )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Subd_Try_PE                                                            */
  /*                                                                         */
  /*  A stripped down version of Subd_Try used in the planar embedding.      */
  /*                                                                         */
  /***************************************************************************/

  int curre;
  int othere;
  int i;
  int j;
  int k;
  int tt;
  int r;

  bool consider_incidente = true;

  while ( consider_incidente )
  {
    consider_incidente = false;

    for ( j = 0; j < M; j++ )
    {

      if ( ( *not_added_edges )[j] == 0 )
      {
        curre = M;
        for ( r = 0; r <= 1; r++ )
        {
          i = ( *Sixidx )[r * M + j] - 1;

          if ( ( *Cidx )[curre] - ( *Cidx )[curre - 1] > 0 )
          {
            if ( ( *Sixidx )[curre - 1] == i + 1 )
            {
              if ( Clabel[( *Cidx )[curre - 1]] == ( *Sixidx )[( 2 * M ) - 1 + curre] )
              {
                othere = ( *Sixidx )[( 2 * M ) - 1 + curre];
                if ( ( *Sixidx )[othere - 1] == i + 1 && Clabel[( *Cidx )[othere - 1]] == curre )
                {
                  consider_incidente = true;
                  for ( k = ( *Cidx )[curre - 1]; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = curre + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }
                  for ( k = ( *Cidx )[othere - 1]; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = othere + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }

                  if ( ( *deg )[i] > 2 )
                  {
                    if ( ( *Sixidx )[( *Sixidx )[( 2 * M ) - 1 + othere] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                    }
                    else
                    {
                      ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                    }
                    if ( ( *Sixidx )[( *Sixidx )[( 3 * M ) - 1 + curre] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                    }
                    else
                    {
                      ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                    }
                    tt = ( *Sixidx )[( 3 * M ) - 1 + curre];
                    ( *Sixidx )[( 2 * M ) - 1 + curre] = ( *Sixidx )[( 2 * M ) - 1 + othere];
                    ( *Sixidx )[( 3 * M ) - 1 + curre] = othere;
                    ( *Sixidx )[( 2 * M ) - 1 + othere] = curre;
                    ( *Sixidx )[( 3 * M ) - 1 + othere] = tt;
                  }
                  current_crossings = current_crossings - 1;
                }
                else if ( ( *Sixidx )[M - 1 + othere] == i + 1 && Clabel[( *Cidx )[othere] - 1] == curre )
                {
                  consider_incidente = true;
                  for ( k = ( *Cidx )[curre - 1]; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = curre + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }
                  for ( k = ( *Cidx )[othere] - 1; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = othere + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }

                  if ( ( *deg )[i] > 2 )
                  {
                    if ( ( *Sixidx )[( *Sixidx )[( 4 * M ) - 1 + othere] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 3 * M ) + ( *Sixidx )[( 4 * M ) - 1 + othere] - 1] = curre;
                    }
                    else
                    {
                      ( *Sixidx )[( 5 * M ) + ( *Sixidx )[( 4 * M ) - 1 + othere] - 1] = curre;
                    }
                    if ( ( *Sixidx )[( *Sixidx )[( 3 * M ) - 1 + curre] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                    }
                    else
                    {
                      ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + curre]] = othere;
                    }
                    tt = ( *Sixidx )[( 3 * M ) - 1 + curre];
                    ( *Sixidx )[( 2 * M ) - 1 + curre] = ( *Sixidx )[( 4 * M ) - 1 + othere];
                    ( *Sixidx )[( 3 * M ) - 1 + curre] = othere;
                    ( *Sixidx )[( 4 * M ) - 1 + othere] = curre;
                    ( *Sixidx )[( 5 * M ) - 1 + othere] = tt;
                  }
                  current_crossings = current_crossings - 1;
                }
              }
            }
            else
            {

              if ( Clabel[( *Cidx )[curre] - 1] == ( *Sixidx )[( 4 * M ) - 1 + curre] )
              {
                othere = ( *Sixidx )[( 4 * M ) - 1 + curre];
                if ( ( *Sixidx )[othere - 1] == i + 1 && Clabel[( *Cidx )[othere - 1]] == curre )
                {
                  consider_incidente = true;
                  for ( k = ( *Cidx )[curre] - 1; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = curre + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }
                  for ( k = ( *Cidx )[othere - 1]; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = othere + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }

                  if ( ( *deg )[i] > 2 )
                  {
                    if ( ( *Sixidx )[( *Sixidx )[( 2 * M ) - 1 + othere] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                    }
                    else
                    {
                      ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + othere]] = curre;
                    }
                    if ( ( *Sixidx )[( *Sixidx )[( 5 * M ) - 1 + curre] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                    }
                    else
                    {
                      ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                    }
                    tt = ( *Sixidx )[( 5 * M ) - 1 + curre];
                    ( *Sixidx )[( 4 * M ) - 1 + curre] = ( *Sixidx )[( 2 * M ) - 1 + othere];
                    ( *Sixidx )[( 5 * M ) - 1 + curre] = othere;
                    ( *Sixidx )[( 2 * M ) - 1 + othere] = curre;
                    ( *Sixidx )[( 3 * M ) - 1 + othere] = tt;
                  }
                  current_crossings = current_crossings - 1;
                }
                else if ( ( *Sixidx )[M - 1 + othere] == i + 1 && Clabel[( *Cidx )[othere] - 1] == curre )
                {
                  consider_incidente = true;
                  for ( k = ( *Cidx )[curre] - 1; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = curre + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }
                  for ( k = ( *Cidx )[othere] - 1; k < ( *Cidx )[M] - 1; k++ )
                  {
                    Clabel[k] = Clabel[k + 1];
                    Dlabel[k] = Dlabel[k + 1];
                  }
                  for ( k = othere + 1; k < M + 2; k++ )
                  {
                    ( *Cidx )[k - 1] = ( *Cidx )[k - 1] - 1;
                  }

                  if ( ( *deg )[i] > 2 )
                  {
                    if ( ( *Sixidx )[( *Sixidx )[( 4 * M ) - 1 + othere] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + othere]] = curre;
                    }
                    else
                    {
                      ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + othere]] = curre;
                    }
                    if ( ( *Sixidx )[( *Sixidx )[( 5 * M ) - 1 + curre] - 1] == i + 1 )
                    {
                      ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                    }
                    else
                    {
                      ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + curre]] = othere;
                    }
                    tt = ( *Sixidx )[( 5 * M ) - 1 + curre];
                    ( *Sixidx )[( 4 * M ) - 1 + curre] = ( *Sixidx )[( 4 * M ) - 1 + othere];
                    ( *Sixidx )[( 5 * M ) - 1 + curre] = othere;
                    ( *Sixidx )[( 4 * M ) - 1 + othere] = curre;
                    ( *Sixidx )[( 5 * M ) - 1 + othere] = tt;
                  }
                  current_crossings = current_crossings - 1;
                }
              }
            }
          }
        }
      }
    }
  }

  int restartloop = 1;
  int* pcnt = (int*)malloc( d * sizeof( int ) );
  int oldN = N;
  int oldM = M;
  int ci;
  int qi;
  int qi1;
  int cj;
  int turnchk1;
  int chksum = 0;
  int* cc = (int*)malloc( maxpdist * sizeof( int ) );
  int v2;
  int* newc1 = (int*)malloc( maxpdist * sizeof( int ) );
  int* newc2 = (int*)malloc( maxpdist * sizeof( int ) );
  int* newd1 = (int*)malloc( maxpdist * sizeof( int ) );
  int* newd2 = (int*)malloc( maxpdist * sizeof( int ) );
  int ncnt1;
  int ncnt2;
  int des;
  int cnt3;
  int cnt4;
  int qc;
  int qf;
  int* pe = (int*)malloc( maxpdist * sizeof( int ) );
  int* ze = (int*)malloc( maxpdist * sizeof( int ) );
  int* de = (int*)malloc( maxpdist * sizeof( int ) );
  int erg;
  int qrg;
  int cntr;
  int cntr2;
  int* qe = (int*)malloc( maxpdist * sizeof( int ) );
  int q;
  int* Mid = (int*)malloc( d * maxpdist * sizeof( int ) );
  int* oe = (int*)malloc( maxpdist * sizeof( int ) );

  while ( restartloop == 1 )
  {
    restartloop = 0;

    for ( i = 1; i < d + M - oldM + 1; i++ )
    {
      if ( i <= d )
      {
        ci = fi[i - 1];
      }
      else
      {
        ci = oldM + i - d;
      }
      qi = ( *Cidx )[ci - 1];
      qi1 = ( *Cidx )[ci];
      if ( ( *Sixidx )[M + ci - 1] > oldN )
      {
        turnchk1 = 1;
      }
      else
      {
        turnchk1 = 0;
      }
      for ( j = 0; j < qi1 - qi; j++ )
      {
        if ( turnchk1 == 1 )
        {
          cj = qi1 - j - 1;
        }
        else
        {
          cj = qi + j;
        }
        chksum = 0;
        for ( k = qi; k < qi1; k++ )
        {
          if ( Clabel[k] == Clabel[cj] )
          {
            chksum = chksum + 1;
          }
        }
        if ( chksum > 1 )
        {
          // then subdvision incoming so reallocate (*Sixidx), (*Cidx), (*incidente), (*deg)
          // fprintf(stdout,"SUBDIVIDING");
          *Sixidx = (int*)realloc( *Sixidx, ( 6 * M + 6 ) * sizeof( int ) );
          for ( r = 0; r < 5; r++ )
          {
            for ( k = r * M; k < ( r + 1 ) * M; k++ )
            {

              ( *Sixidx )[6 * M + 6 - k - 2 - r] = ( *Sixidx )[6 * M - k - 1];
            }
          }
          for ( k = 0; k < 5; k++ )
          {
            ( *Sixidx )[( k + 1 ) * ( M + 1 ) - 1] = 0;
          }
          *incidente = (int*)realloc( *incidente, ( N + 1 ) * sizeof( int ) );
          *deg = (int*)realloc( *deg, ( N + 1 ) * sizeof( int ) );
          *Cidx = (int*)realloc( *Cidx, ( M + 2 ) * sizeof( int ) );
          *checked = (int*)realloc( *checked, ( N + 1 ) * sizeof( int ) );
          *not_added_edges = (int*)realloc( *not_added_edges, ( M + 1 ) * sizeof( int ) );

          for ( k = 0; k < maxpdist; k++ )
          { // reset cc for use again
            cc[k] = 0;
          }
          v2 = Clabel[cj];
          N = N + 1;
          M = M + 1;
          cntr = 0;
          cntr2 = 0;
          ncnt1 = -1;
          ncnt2 = -1;

          if ( ( *Sixidx )[M - 1 + ci] <= oldN )
          {
            for ( k = qi; k < cj + 1; k++ )
            {
              ncnt1 = ncnt1 + 1;
              newc1[ncnt1] = Clabel[k];
              newd1[ncnt1] = Dlabel[k];
            }
            for ( k = 0; k < ( qi1 - cj - 1 ); k++ )
            {

              newc2[k] = Clabel[qi1 - 1 - k];
              newd2[k] = -Dlabel[qi1 - 1 - k];
            }
            ncnt2 = ( qi1 - cj - 2 );
          }
          else
          {
            for ( k = 0; k < ( qi1 - cj ); k++ )
            {
              ncnt1 = ncnt1 + 1;
              newc1[k] = Clabel[qi1 - 1 - k];
              newd1[k] = -Dlabel[qi1 - 1 - k];
            }
            for ( k = qi; k < cj; k++ )
            {

              ncnt2 = ncnt2 + 1;
              newc2[ncnt2] = Clabel[k];
              newd2[ncnt2] = Dlabel[k];
            }
          }
          cntr = -1;
          qc = ( *Cidx )[v2 - 1];
          for ( k = qc; k < ( *Cidx )[v2]; k++ )
          {
            if ( Clabel[k] == ci )
            {
              cntr = cntr + 1;
              cc[cntr] = k - qc;
            }
          }

          for ( k = 0; k < maxpdist; k++ )
          {
            pe[k] = 0;
            qe[k] = 0;
            ze[k] = 0;
          }

          cntr = -1;
          if ( ( *Sixidx )[M - 1 + ci] > oldN )
          {
            q = Mid[ci - oldM - 1];
            if ( ( *Sixidx )[ci - 1] == fin[q - 1] )
            {
              for ( k = pcnt[q - 1] + 1; k < pdist[q - 1] + 1; k++ )
              {
                cntr = cntr + 1;
                pe[cntr] = pathedge[pidx[q - 1] + k - 1];
              }
            }
            else
            {
              for ( k = 1; k < pcnt[q - 1]; k++ )
              {
                cntr = cntr + 1;
                pe[cntr] = pathedge[pidx[q - 1] + pcnt[q - 1] - k - 1];
              }
            }
          }
          else
          {
            for ( k = 0; k < d; k++ )
            {
              if ( fi[k] == ci )
              {
                qf = k;
                break;
              }
            }
            if ( ( *Sixidx )[ci - 1] == vm )
            {
              for ( k = pidx[qf]; k < pidx[qf + 1]; k++ )
              {
                cntr = cntr + 1;
                pe[cntr] = pathedge[k];
              }
            }
            else
            {
              for ( k = 0; k < pdist[qf]; k++ )
              {
                cntr = cntr + 1;
                pe[cntr] = pathedge[pidx[qf + 1] - 1 - k];
              }
            }
          }

          cntr2 = -1;
          for ( k = 0; k < cntr + 1; k++ )
          {
            oe[k] = origedgelabel[pe[k] - 1];
            if ( oe[k] == v2 )
            {

              cntr2 = cntr2 + 1;
              ze[cntr2] = k + 1;
              de[cntr2] = pe[k];
            }
          }

          des = de[0];

          qsort( de, cntr2 + 1, sizeof( int ), ascend_cmpfunc );

          if ( ( *Sixidx )[M - 1 + ci] > oldN )
          {
            Mid[M - oldM - 1] = q;
            if ( ( *Sixidx )[ci - 1] == fin[q - 1] )
            {
              pcnt[q - 1] = pcnt[q - 1] + ze[0];
            }
            else
            {
              pcnt[q - 1] = pcnt[q - 1] - ze[0];
            }
          }
          else
          {

            Mid[M - oldM - 1] = qf + 1;
            if ( ( *Sixidx )[ci - 1] == vm )
            {
              pcnt[qf] = ze[0];
            }
            else
            {
              pcnt[qf] = pdist[qf] + 1 - ze[0];
            }
          }

          bool lessN = false;
          if ( ( *Sixidx )[M - 1 + ci] <= oldN )
            lessN = true;

          for ( k = 0; k < cntr2 + 1; k++ )
          {
            if ( de[k] != des )
            {
              Clabel[qc + cc[k]] = M;
              if ( lessN )
              {
                Dlabel[qc + cc[k]] = -Dlabel[qc + cc[k]];
              }
            }
            else
            {
              if ( lessN == false )
              {
                Dlabel[qc + cc[k]] = -Dlabel[qc + cc[k]];
              }
            }
          }

          if ( ( *Sixidx )[M - 1 + ci] <= oldN )
          {
            ( *Sixidx )[M - 1] = ( *Sixidx )[M - 1 + ci];
            ( *Sixidx )[( 2 * M ) - 1] = N;
            ( *Sixidx )[( 3 * M ) - 1] = ( *Sixidx )[( 4 * M ) - 1 + ci];
            ( *Sixidx )[( 4 * M ) - 1] = ( *Sixidx )[( 5 * M ) - 1 + ci];
            ( *Sixidx )[( 5 * M ) - 1] = ci;
            ( *Sixidx )[( 6 * M ) - 1] = ci;
            if ( ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + ci]] == ci )
            {
              ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + ci]] = M;
            }
            else
            {
              ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 4 * M ) - 1 + ci]] = M;
            }
            if ( ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + ci]] == ci )
            {
              ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + ci]] = M;
            }
            else
            {
              ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 5 * M ) - 1 + ci]] = M;
            }
            ( *Sixidx )[M - 1 + ci] = N;
            ( *Sixidx )[( 4 * M ) - 1 + ci] = M;
            ( *Sixidx )[( 5 * M ) - 1 + ci] = M;
          }
          else
          {
            ( *Sixidx )[M - 1] = ( *Sixidx )[ci - 1];
            ( *Sixidx )[( 2 * M ) - 1] = N;
            ( *Sixidx )[( 3 * M ) - 1] = ( *Sixidx )[( 2 * M ) - 1 + ci];
            ( *Sixidx )[( 4 * M ) - 1] = ( *Sixidx )[( 3 * M ) - 1 + ci];
            ( *Sixidx )[( 5 * M ) - 1] = ci;
            ( *Sixidx )[( 6 * M ) - 1] = ci;
            if ( ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + ci]] == ci )
            {
              ( *Sixidx )[( 3 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + ci]] = M;
            }
            else
            {
              ( *Sixidx )[( 5 * M ) - 1 + ( *Sixidx )[( 2 * M ) - 1 + ci]] = M;
            }
            if ( ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + ci]] == ci )
            {
              ( *Sixidx )[( 2 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + ci]] = M;
            }
            else
            {
              ( *Sixidx )[( 4 * M ) - 1 + ( *Sixidx )[( 3 * M ) - 1 + ci]] = M;
            }
            ( *Sixidx )[ci - 1] = ( *Sixidx )[M - 1 + ci];
            ( *Sixidx )[M - 1 + ci] = N;
            ( *Sixidx )[( 2 * M ) - 1 + ci] = ( *Sixidx )[( 4 * M ) - 1 + ci];
            ( *Sixidx )[( 3 * M ) - 1 + ci] = ( *Sixidx )[( 4 * M ) - 1 + ci];
            ( *Sixidx )[( 4 * M ) - 1 + ci] = M;
            ( *Sixidx )[( 5 * M ) - 1 + ci] = M;
          }

          for ( k = 0; k < ncnt2 + 1; k++ )
          {
            if ( newc2[k] != v2 )
            {
              erg = newc2[k];
              for ( r = ( *Cidx )[erg - 1]; r < ( *Cidx )[erg]; r++ )
              {
                if ( Clabel[r] == ci )
                {

                  Clabel[r] = M;
                  if ( lessN )
                    Dlabel[r] = -Dlabel[r];
                }
              }
            }
          }
          if ( lessN == false )
          {
            for ( k = 0; k < ncnt1 + 1; k++ )
            {
              if ( newc1[k] != v2 )
              {
                erg = newc1[k];
                for ( r = ( *Cidx )[erg - 1]; r < ( *Cidx )[erg]; r++ )
                {
                  if ( Clabel[r] == ci )
                  {

                    Dlabel[r] = -Dlabel[r];
                  }
                }
              }
            }
          }
          cnt3 = -1;
          cnt4 = -1;
          for ( k = 0; k < ( *Cidx )[M - 1]; k++ )
          {
            if ( k >= qi && k <= qi + ncnt1 )
            { // these need checked
              cnt3 = cnt3 + 1;
              Clabel[k] = newc1[cnt3];
              Dlabel[k] = newd1[cnt3];
            }

            if ( k > qi + ncnt1 && k < ( *Cidx )[M - 1] - ncnt2 - 1 )
            {
              Clabel[k] = Clabel[k + ncnt2 + 1];
              Dlabel[k] = Dlabel[k + ncnt2 + 1];
            }
            if ( k >= ( *Cidx )[M - 1] - ncnt2 - 1 )
            {
              cnt4 = cnt4 + 1;
              Clabel[k] = newc2[cnt4];
              Dlabel[k] = newd2[cnt4];
            }
          }
          for ( k = ci; k < M; k++ )
          {
            ( *Cidx )[k] = ( *Cidx )[k] - cnt4 - 1;
          }
          ( *Cidx )[M] = ( *Cidx )[M - 1] + cnt4 + 1;
          ( *incidente )[N - 1] = M;
          ( *incidente )[( *Sixidx )[M - 1] - 1] = M;
          ( *deg )[N - 1] = 2;
          ( *checked )[N - 1] = 1;
          ( *not_added_edges )[M - 1] = 0;

          restartloop = 1;

          break;
        }
      }
      if ( restartloop == 1 )
      {
        break;
      }
    }
  }

  outputs[0] = N;
  outputs[1] = M;
  outputs[2] = current_crossings;

  free( pcnt );
  free( cc );
  free( newc1 );
  free( newc2 );
  free( newd1 );
  free( newd2 );
  free( pe );
  free( ze );
  free( de );
  free( qe );
  free( Mid );
  free( oe );
}

void New_Clabel_PE( int M, int Md, int* Sixidx_plan, int* Sixidx, int* Clabel, int* Dlabel, int* Cidx, int* elist, int* eidx, int* pathedge, int* paths, int* pidx, int* pathsidx, int* pdist, int maxpdist, int d, int* fin, int* fi, int* F, int* Fdir, int* Fidx, int vm, int* origedgelabel, int* origedgeorder, int newface, int* pathedgeonce, int treecheck, int maxdeadjm )
{

  /***************************************************************************/
  /*                                                                         */
  /*  New_Clabel_PE                                                          */
  /*                                                                         */
  /*  A stripped down version of New_Clabel used in the planar embedding.    */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int change;
  int comb_ancestor_count = 0;
  int* comb_ancestor_sparse;
  comb_ancestor_sparse = (int*)malloc( Md * sizeof( int ) );
  int* comb_ancestor_first;
  comb_ancestor_first = (int*)malloc( Md * sizeof( int ) );
  int* comb_ancestor_next;
  comb_ancestor_next = (int*)malloc( Md * sizeof( int ) );
  int* eparent;
  eparent = (int*)malloc( Md * sizeof( int ) );
  int* vparent;
  vparent = (int*)malloc( Md * sizeof( int ) );
  int* comb_count;
  comb_count = (int*)malloc( Md * sizeof( int ) );
  int* comb_ancestor_recent;
  comb_ancestor_recent = (int*)malloc( Md * sizeof( int ) );

  if ( treecheck == 1 )
  {
    for ( i = 0; i < Md; i++ )
    {
      comb_ancestor_sparse[i] = 0;
      comb_ancestor_first[i] = 0;
      comb_ancestor_next[i] = 0;
      eparent[i] = 0;
      vparent[i] = 0;
      comb_count[i] = 0;
      comb_ancestor_recent[i] = 0;
    }
    comb_ancestor_count = Paths_Tree( Md, pidx, pathedge, paths, Sixidx_plan, d, fin, F, Fdir, Fidx, pdist, maxpdist, pathsidx, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, eparent, vparent, comb_count, comb_ancestor_recent );
  }

  int sumdist = 0;
  for ( i = 0; i < d; i++ )
  {
    sumdist = sumdist + pdist[i];
  }

  int* newClabel;
  int* newDlabel;
  newClabel = (int*)malloc( ( 2 * sumdist + Cidx[M] ) * sizeof( int ) );
  newDlabel = (int*)malloc( ( 2 * sumdist + Cidx[M] ) * sizeof( int ) );
  for ( i = 0; i < Cidx[M]; i++ )
  {
    newClabel[i] = 0;
    newDlabel[i] = 0;
  }

  int* newCidx = (int*)malloc( ( M + 1 ) * sizeof( int ) );
  for ( i = 0; i < M + 1; i++ )
  {
    newCidx[i] = 0;
  }
  int ci;
  int q;
  int qn;
  int k;
  int r;
  int x1cnt = -1;
  int y1cnt = -1;
  int* x1;
  x1 = (int*)malloc( d * sizeof( int ) );
  int* y1;
  y1 = (int*)malloc( d * sizeof( int ) );
  int same;
  int backe;
  int epos;
  int f1;
  int last_i = 0;
  int ordcnt;
  int* order = (int*)malloc( d * sizeof( int ) );
  int face;
  int* ctemp;
  ctemp = (int*)malloc( M * sizeof( int ) ); // max crossings size, but M for now
  for ( i = 0; i < M; i++ )
    ctemp[i] = 0;
  int ctempcnt = -1;
  int ori;

  int dtempcnt = -1;
  int* dtemp;
  dtemp = (int*)malloc( M * sizeof( int ) );

  int* fi_check;
  fi_check = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
  {
    fi_check[i] = 0;
  }
  for ( i = 0; i < d; i++ )
  {
    fi_check[fi[i] - 1] = 1;
  }

  int newccnt = -1;
  int ec;
  int curre;
  for ( i = 0; i < M; i++ )
  {
    ci = 1;
    for ( j = 0; j < d; j++ )
    {
      if ( i + 1 == fi[j] )
      {
        ci = 0;
        break;
      }
    }

    if ( ci == 1 )
    {
      ctempcnt = -1;
      dtempcnt = -1;
      x1cnt = -1;
      y1cnt = -1;
      qn = eidx[i + 1] - eidx[i];
      for ( j = 0; j < qn; j++ )
      {
        q = elist[eidx[i] + j];
        if ( treecheck == 1 && comb_count[q - 1] )
        {

          ordcnt = Traverse_Paths_Tree( Md, eparent, vparent, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, q, d, fin, pidx, pathedge, fi, order );
          x1cnt = -1;
          y1cnt = -1;
          for ( k = 0; k < d; k++ )
          {
            for ( r = pidx[k]; r < pidx[k + 1]; r++ )
            {
              if ( pathedge[r] == q )
              {
                x1cnt = x1cnt + 1;
                y1cnt = y1cnt + 1;
                x1[x1cnt] = k + 1;
                y1[y1cnt] = r - pidx[k] + 1;
              }
            }
          }

          same = 1;
          if ( origedgeorder[q - 1] == 1 || origedgeorder[q - 1] == ( eidx[origedgelabel[q - 1]] - eidx[origedgelabel[q - 1] - 1] ) )
          {
            if ( Sixidx[origedgelabel[q - 1] - 1] != Sixidx_plan[q - 1] && Sixidx[M + origedgelabel[q - 1] - 1] != Sixidx_plan[Md + q - 1] )
            {
              same = 0;
            }
          }
          else
          {
            backe = elist[eidx[origedgelabel[q - 1] - 1] + origedgeorder[q - 1] - 2];
            if ( Sixidx_plan[backe - 1] != Sixidx_plan[q - 1] && Sixidx_plan[Md + backe - 1] != Sixidx_plan[q - 1] )
            {
              same = 0;
            }
          }
          for ( k = Fidx[paths[pathsidx[x1[0] - 1] + y1[0]] - 1]; k < Fidx[paths[pathsidx[x1[0] - 1] + y1[0]]]; k++ )
          {
            if ( F[k] == q )
            {
              epos = k;
              break;
            }
          }
          if ( Fdir[epos] == 1 )
          {
            if ( same == 1 )
            {

              same = 0;
            }
            else
            {

              same = 1;
            }
          }

          if ( same == 0 )
          {
            for ( k = 0; k < ordcnt; k++ )
            {
              ctempcnt = ctempcnt + 1;
              ctemp[ctempcnt] = order[ordcnt - k - 1];
            }
            for ( k = 0; k < ordcnt; k++ )
            {
              for ( r = 0; r < d; r++ )
              {
                if ( fi[r] == order[ordcnt - k - 1] )
                {
                  f1 = r + 1;
                  break;
                }
              }
              face = paths[pathsidx[x1[x1cnt - k] - 1] + y1[y1cnt - k]];
              ori = New_Crossing_Orientation( Md, Sixidx_plan, Fidx[face] - Fidx[face - 1], F, Fidx, origedgeorder, q, vm, fin[f1 - 1], face, origedgelabel[q - 1], elist, eidx );

              dtempcnt = dtempcnt + 1;
              dtemp[dtempcnt] = ori;
            }
          }
          else
          {

            for ( k = 0; k < ordcnt; k++ )
            {
              ctempcnt = ctempcnt + 1;

              ctemp[ctempcnt] = order[k];
            }
            for ( k = 0; k < ordcnt; k++ )
            {
              for ( r = 0; r < d; r++ )
              {
                if ( fi[r] == order[k] )
                {
                  f1 = r + 1;
                  break;
                }
              }
              face = paths[pathsidx[x1[k] - 1] + y1[k]];
              ori = New_Crossing_Orientation( Md, Sixidx_plan, Fidx[face] - Fidx[face - 1], F, Fidx, origedgeorder, q, vm, fin[f1 - 1], face, origedgelabel[q - 1], elist, eidx );
              dtempcnt = dtempcnt + 1;
              dtemp[dtempcnt] = ori;
            }
          }
        }
        else
        {
          x1[0] = pathedgeonce[q - 1];
          y1[0] = pathedgeonce[maxdeadjm + q - 1]; // check this length of pathedgeonce.
          if ( x1[0] != 0 )
          {
            ctempcnt = ctempcnt + 1;
            ctemp[ctempcnt] = fi[x1[0] - 1];

            face = paths[pathsidx[x1[0] - 1] + y1[0]];
            ori = New_Crossing_Orientation( Md, Sixidx_plan, Fidx[face] - Fidx[face - 1], F, Fidx, origedgeorder, q, vm, fin[x1[0] - 1], face, origedgelabel[q - 1], elist, eidx );
            dtempcnt = dtempcnt + 1;
            dtemp[dtempcnt] = ori;
          }
        }

        if ( j < qn - 1 && fi_check[Clabel[Cidx[i] + j] - 1] == 0 )
        {
          ctempcnt = ctempcnt + 1;
          ctemp[ctempcnt] = Clabel[Cidx[i] + j];
          dtempcnt = dtempcnt + 1;
          dtemp[dtempcnt] = Dlabel[Cidx[i] + j];
        }
      }

      if ( ctempcnt != -1 )
      {
        for ( k = 0; k < ctempcnt + 1; k++ )
        {
          newccnt = newccnt + 1;
          newClabel[newccnt] = ctemp[k];
          newDlabel[newccnt] = dtemp[k];
        }
      }

      if ( last_i != 0 )
      {
        for ( k = last_i + 2; k <= i + 1; k++ )
        {
          newCidx[k - 1] = newCidx[last_i];
        }
      }
      last_i = i + 1;
      newCidx[i + 1] = newCidx[i] + ctempcnt + 1;

      if ( i != M - 1 )
      {
        for ( k = i + 3; k <= M + 1; k++ )
        {
          newCidx[k - 1] = newCidx[i + 1];
        }
      }
    }
  }
  free( ctemp );
  free( dtemp );
  free( fi_check );

  for ( i = 0; i < newCidx[M]; i++ )
  {

    Clabel[i] = newClabel[i];
    Dlabel[i] = newDlabel[i];
  }

  free( newClabel );
  free( newDlabel );
  for ( i = 0; i < M + 1; i++ )
  {
    Cidx[i] = newCidx[i];
  }

  free( newCidx );

  int nf = Fidx[newface] - Fidx[newface - 1];
  int* fincheck = (int*)malloc( d * sizeof( int ) );
  for ( i = 0; i < d; i++ )
  {
    fincheck[i] = 0;
  }

  int* atemp = (int*)malloc( d * sizeof( int ) );
  int vstart;
  int acnt = -1;
  int nexteidx;
  x1cnt = -1;
  for ( i = 0; i < nf; i++ )
  {
    q = -1;
    if ( Fdir[Fidx[newface - 1] + i] == 1 )
    {
      vstart = Sixidx_plan[F[Fidx[newface - 1] + i] - 1];
    }
    else
    {
      vstart = Sixidx_plan[Md + F[Fidx[newface - 1] + i] - 1];
    }
    for ( j = 0; j < d; j++ )
    {
      if ( fin[j] == vstart )
      {
        q = j;
        break;
      }
    }
    if ( q != -1 && fincheck[q] == 0 )
    {

      acnt = acnt + 1;
      atemp[acnt] = fi[q];
      fincheck[q] = 1;
      curre = fi[q];
      if ( i > 0 )
      {
        nexteidx = Fidx[newface - 1] + i - 1;
      }
      else
      {
        nexteidx = Fidx[newface] - 1;
      }
      if ( Sixidx[origedgelabel[curre - 1] - 1] == fin[q] )
      {
        Sixidx[( 2 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[nexteidx] - 1];
        Sixidx[( 3 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[Fidx[newface - 1] + i] - 1];
      }
      else
      {
        Sixidx[( 4 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[nexteidx] - 1];
        Sixidx[( 5 * M ) + origedgelabel[curre - 1] - 1] = origedgelabel[F[Fidx[newface - 1] + i] - 1];
      }
      if ( Sixidx[origedgelabel[F[Fidx[newface - 1] + i] - 1] - 1] == fin[q] )
      {
        Sixidx[( 2 * M ) + origedgelabel[F[Fidx[newface - 1] + i] - 1] - 1] = fi[q];
      }
      else
      {
        Sixidx[( 4 * M ) + origedgelabel[F[Fidx[newface - 1] + i] - 1] - 1] = fi[q];
      }
      if ( Sixidx[origedgelabel[F[nexteidx] - 1] - 1] == fin[q] )
      {
        Sixidx[( 3 * M ) + origedgelabel[F[nexteidx] - 1] - 1] = fi[q];
      }
      else
      {
        Sixidx[( 5 * M ) + origedgelabel[F[nexteidx] - 1] - 1] = fi[q];
      }
    }
    x1cnt = -1;
    for ( j = 0; j < d; j++ )
    {
      if ( pdist[j] > 0 && pathedge[pidx[j]] == F[Fidx[newface - 1] + i] )
      {
        x1cnt = x1cnt + 1;
        x1[x1cnt] = j + 1;
      }
    }

    if ( x1cnt > 0 )
    {
      ordcnt = Traverse_Paths_Tree( Md, eparent, vparent, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, F[Fidx[newface - 1] + i], d, fin, pidx, pathedge, fi, order );
      for ( j = 0; j < ordcnt; j++ )
      {
        acnt = acnt + 1;
        atemp[acnt] = order[j];
      }
    }
    else if ( x1cnt == 0 )
    {
      acnt = acnt + 1;
      atemp[acnt] = fi[x1[0] - 1];
    }
  }

  free( fincheck );

  int fs;
  int bs;
  for ( i = 0; i < d; i++ )
  {
    if ( i == 0 )
    {
      bs = atemp[acnt];
    }
    else
    {
      bs = atemp[i - 1];
    }
    if ( i == d - 1 )
    {
      fs = atemp[0];
    }
    else
    {
      fs = atemp[i + 1];
    }
    if ( Sixidx[atemp[i] - 1] == vm )
    {
      Sixidx[( 2 * M ) + atemp[i] - 1] = fs;
      Sixidx[( 3 * M ) + atemp[i] - 1] = bs;
    }
    else
    {
      Sixidx[( 4 * M ) + atemp[i] - 1] = fs;
      Sixidx[( 5 * M ) + atemp[i] - 1] = bs;
    }
  }
  free( atemp );
  int* x12 = (int*)malloc( maxpdist * sizeof( int ) );
  int* checked;
  int* newedgelabel;
  int* newedgelabelun;
  if ( maxpdist > 0 )
  {
    checked = (int*)malloc( maxpdist * sizeof( int ) );
    newedgelabel = (int*)malloc( maxpdist * sizeof( int ) );
    newedgelabelun = (int*)malloc( maxpdist * sizeof( int ) );
  }
  int tj;
  int sr;

  for ( i = 0; i < d; i++ )
  {
    fs = -1;
    if ( pdist[i] > 0 )
    {
      x1cnt = -1;
      fs = -1;
      for ( j = Cidx[fi[i] - 1]; j < Cidx[M]; j++ )
      {
        fs = fs + 1;
        Clabel[Cidx[M] - 1 + pdist[i] - fs] = Clabel[Cidx[M] - 1 - fs];
      }
      fs = -1;

      for ( j = Cidx[fi[i] - 1]; j < Cidx[M]; j++ )
      {
        fs = fs + 1;

        Dlabel[Cidx[M] - 1 + pdist[i] - fs] = Dlabel[Cidx[M] - 1 - fs];
      }
      fs = -1;
      if ( vm < fin[i] )
      {
        for ( j = Cidx[fi[i] - 1]; j < Cidx[fi[i] - 1] + pdist[i]; j++ )
        {
          fs = fs + 1;
          Clabel[j] = origedgelabel[pathedge[pidx[i] + fs] - 1];
        }
      }
      else
      {
        for ( j = Cidx[fi[i] - 1]; j < Cidx[fi[i] - 1] + pdist[i]; j++ )
        {
          fs = fs + 1;
          Clabel[j] = origedgelabel[pathedge[pidx[i + 1] - 1 - fs] - 1];
        }
      }
      for ( j = fi[i]; j <= M; j++ )
      {
        Cidx[j] = Cidx[j] + pdist[i];
      }
      for ( j = 0; j < maxpdist; j++ )
      {
        checked[j] = 0;
      }
      for ( j = 0; j < pdist[i]; j++ )
      {
        x1cnt = -1;
        for ( k = 0; k < pdist[i]; k++ )
        {
          if ( Clabel[Cidx[fi[i] - 1] + j] == Clabel[Cidx[fi[i] - 1] + k] )
          {
            x1cnt = x1cnt + 1;
            x12[x1cnt] = k + 1;
          }
        }
        if ( x1cnt > 0 && checked[j] == 0 )
        {

          for ( k = 0; k <= x1cnt; k++ )
          {
            checked[x12[k] - 1] = 1;
            if ( vm < fin[i] )
            {
              newedgelabel[k] = pathedge[pidx[i] + x12[k] - 1];
            }
            else
            {
              newedgelabel[k] = pathedge[pidx[i + 1] - x12[k]];
            }
            newedgelabelun[k] = newedgelabel[k];
          }
          qsort( newedgelabel, x1cnt + 1, sizeof( int ), ascend_cmpfunc );
          tj = Clabel[Cidx[fi[i] - 1] + j];
          fs = -1;
          for ( k = Cidx[tj - 1]; k < Cidx[tj]; k++ )
          {
            if ( Clabel[k] == fi[i] )
            {
              fs = fs + 1;
              sr = -1;
              for ( r = 0; r < x1cnt + 1; r++ )
              {
                if ( newedgelabelun[r] == newedgelabel[fs] )
                {
                  sr = r;
                }
              }
              Dlabel[Cidx[fi[i] - 1] + x12[sr] - 1] = -Dlabel[k]; // may be wrong here.
            }
          }
        }
        else if ( x1cnt == 0 )
        {
          tj = Clabel[Cidx[fi[i] - 1] + j];
          for ( k = Cidx[tj - 1]; k < Cidx[tj]; k++ )
          {

            if ( Clabel[k] == fi[i] )
            {
              Dlabel[Cidx[fi[i] - 1] + j] = -Dlabel[k];
              break;
            }
          }
        }
      }
    }
  }
  free( x12 );
  free( x1 );
  free( y1 );
  if ( maxpdist > 0 )
  {
    free( newedgelabel );
    free( newedgelabelun );
    free( checked );
  }

  int* checked2 = (int*)malloc( d * sizeof( int ) );
  for ( j = 0; j < d; j++ )
  {
    checked2[j] = 0;
  }
  int* endfaces = (int*)malloc( d * sizeof( int ) );
  int startedge;
  int starte;
  int fn;
  int count;
  int fk;
  int ccnt = -1;
  int* ctemp2 = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
    ctemp2[i] = 150;
  int facee;

  for ( j = 0; j < d; j++ )
    endfaces[j] = paths[pathsidx[j + 1] - 1];

  for ( i = 0; i < d; i++ )
  {
    if ( checked2[i] == 0 )
    {
      for ( j = 0; j < d; j++ )
      {
        if ( endfaces[j] == endfaces[i] )
        {
          checked2[j] = 1;
        }
      }
      if ( pdist[i] > 0 )
      {
        startedge = pathedge[pidx[i + 1] - 1];
        facee = endfaces[i];
        ccnt = -1;

        if ( treecheck > 0 && comb_count[startedge - 1] > 0 )
        {
          ordcnt = Traverse_Paths_Tree( Md, eparent, vparent, comb_ancestor_sparse, comb_ancestor_first, comb_ancestor_next, startedge, d, fin, pidx, pathedge, fi, order );
          for ( k = 0; k < ordcnt; k++ )
          {
            ccnt = ccnt + 1;
            ctemp2[ccnt] = order[k];
          }
        }
        else
        {
          ccnt = ccnt + 1;
          ctemp2[ccnt] = fi[i];
        }
        for ( j = Fidx[facee - 1]; j < Fidx[facee]; j++ )
        {
          if ( F[j] == startedge )
          {
            fs = j;
            break;
          }
        }
        fn = Fidx[facee] - Fidx[facee - 1];
        count = 1;
        for ( j = 0; j < fn; j++ )
        {
          if ( fs < Fidx[facee] - 1 )
          {
            fs = fs + 1;
          }
          else
          {
            fs = Fidx[facee - 1];
          }

          if ( fs > Fidx[facee - 1] )
          {
            bs = fs - 1;
          }
          else
          {
            bs = Fidx[facee] - 1;
          }

          if ( Fdir[fs] == 1 )
          {
            vstart = Sixidx_plan[F[fs] - 1];
          }
          else
          {
            vstart = Sixidx_plan[Md + F[fs] - 1];
          }
          x1cnt = 0;
          fk = -1;
          for ( k = 0; k < d; k++ )
          {
            if ( fin[k] == vstart )
            {
              fk = k;
              x1cnt = x1cnt + 1;
              if ( endfaces[k] != facee )
              {
                x1cnt = -1;
                break;
              }
            }
          }
          if ( fk != -1 && x1cnt > 0 )
          {
            if ( ctemp2[count - 1] == fi[fk] )
            {
              count = count + 1;
              if ( Sixidx[fi[fk] - 1] == vstart )
              {
                Sixidx[( 2 * M ) + fi[fk] - 1] = origedgelabel[F[bs] - 1];
                Sixidx[( 3 * M ) + fi[fk] - 1] = origedgelabel[F[fs] - 1];
              }
              else
              {
                Sixidx[( 4 * M ) + fi[fk] - 1] = origedgelabel[F[bs] - 1];
                Sixidx[( 5 * M ) + fi[fk] - 1] = origedgelabel[F[fs] - 1];
              }
              if ( Sixidx[origedgelabel[F[fs] - 1] - 1] == vstart )
              {
                Sixidx[( 2 * M ) + origedgelabel[F[fs] - 1] - 1] = fi[fk];
              }
              else
              {
                Sixidx[( 4 * M ) + origedgelabel[F[fs] - 1] - 1] = fi[fk];
              }
              if ( Sixidx[origedgelabel[F[bs] - 1] - 1] == vstart )
              {
                Sixidx[( 3 * M ) + origedgelabel[F[bs] - 1] - 1] = fi[fk];
              }
              else
              {
                Sixidx[( 5 * M ) + origedgelabel[F[bs] - 1] - 1] = fi[fk];
              }
            }
          }
          x1cnt = 0;
          for ( k = 0; k < pidx[d]; k++ )
          {
            if ( pathedge[k] == F[fs] )
            {
              x1cnt = x1cnt + 1;
            }
          }
          if ( x1cnt > 0 && F[fs] != startedge )
          {
            count = count + x1cnt;
          }
          if ( count > ccnt + 1 )
          {

            break;
          }
        }
      }
    }
  }

  free( endfaces );
  free( ctemp2 );
  free( comb_ancestor_sparse );
  free( comb_ancestor_first );
  free( comb_ancestor_next );
  free( eparent );
  free( vparent );
  free( comb_count );
  free( checked2 );
  free( comb_ancestor_recent );
  free( order );
}

void Get_Distances_PE( int* dvadjm, int* Fidxm, int Ndm, int Mdm, int M_plan, int* Sixidx, int* dualel, int requiredimprovement, int d, int* fin, int* incidente, int* deg, int current_crossings_best, int* deadjm, int maxdegfin, int maxdeadjm, int** pathedge, int* pidx, int* pathedgeonce, int** paths, int* pathsidx, int* pdist, int* delv, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Get_Distances_PE                                                       */
  /*                                                                         */
  /*  A stripped down version of Get_Distances used in the planar embedding. */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int k;
  int j;
  int newcross = 0;
  int newface = 0;
  int improvementfound = 0;
  int maxpdist = 0;
  int treecheck = 0;
  int used_bigface = 0;

  int Nd = Ndm + 1;

  int* dist = (int*)malloc( d * Nd * sizeof( int ) );
  int* prev = (int*)malloc( d * Nd * sizeof( int ) );
  for ( i = 0; i < d * Nd; i++ )
  {
    dist[i] = 0;
    prev[i] = 0;
  }

  for ( k = 1; k <= d; k++ )
  {
    int* dvadj = (int*)malloc( 3 * Mdm * sizeof( int ) );
    int kkk;
    for ( kkk = 0; kkk < 2 * Mdm; kkk++ )
      dvadj[kkk] = dvadjm[kkk];
    for ( kkk = 2 * Mdm; kkk < 3 * Mdm; kkk++ )
      dvadj[kkk] = 0;

    int* Fidx = (int*)malloc( ( Ndm + 2 ) * sizeof( int ) );
    for ( kkk = 0; kkk < Ndm + 1; kkk++ )
      Fidx[kkk] = Fidxm[kkk];

    int curre = incidente[fin[k - 1] - 1];

    int nexte = -1;
    if ( Sixidx[curre - 1] == fin[k - 1] )
      nexte = Sixidx[curre - 1 + 2 * M_plan];
    else
      nexte = Sixidx[curre - 1 + 4 * M_plan];

    int* newv = (int*)malloc( Mdm * sizeof( int ) );
    for ( i = 0; i < Mdm; i++ )
      newv[i] = 0;
    int cnt = 0;

    int temp_curre = curre;
    int temp_nexte = nexte;

    for ( i = 1; i <= deg[fin[k - 1] - 1]; i++ )
    {
      int fc = -1;
      if ( Sixidx[temp_curre - 1] == fin[k - 1] )
        fc = dualel[temp_curre - 1];
      else
        fc = dualel[temp_curre - 1 + M_plan];

      delv[k - 1 + ( i - 1 ) * d] = fc;

      temp_curre = temp_nexte;

      if ( Sixidx[temp_curre - 1] == fin[k - 1] )
        temp_nexte = Sixidx[temp_curre - 1 + 2 * M_plan];
      else
        temp_nexte = Sixidx[temp_curre - 1 + 4 * M_plan];
    }
    for ( i = 1; i <= deg[fin[k - 1] - 1]; i++ )
    {

      int fc = -1;
      if ( Sixidx[curre - 1] == fin[k - 1] )
        fc = dualel[curre - 1];
      else
        fc = dualel[curre - 1 + M_plan];

      for ( j = 1; j <= Fidx[fc] - Fidx[fc - 1]; j++ )
      {

        bool createEdge = true;

        int jj;
        for ( jj = 1; jj <= maxdegfin; jj++ )
        {
          if ( delv[k - 1 + ( jj - 1 ) * d] == dvadj[Fidx[fc - 1] + j - 1] )
          {
            createEdge = false;
            break;
          }
        }

        if ( createEdge )
        {
          cnt = cnt + 1;
          newv[cnt - 1] = dvadj[Fidx[fc - 1] + j - 1];
          int* f1 = (int*)malloc( ( 8 * ( Fidx[newv[cnt - 1]] - Fidx[newv[cnt - 1] - 1] ) ) * sizeof( int ) );
          int f1cnt = 0;
          int r;
          for ( r = Fidx[newv[cnt - 1] - 1] + 1; r <= Fidx[newv[cnt - 1]]; r++ )
          {
            if ( dvadj[r - 1] == fc )
            {
              f1[f1cnt] = r - Fidx[newv[cnt - 1] - 1];
              f1cnt = f1cnt + 1;
            }
          }

          for ( r = 1; r <= f1cnt; r++ )
          {
            dvadj[Fidx[newv[cnt - 1] - 1] + f1[r - 1] - 1] = Nd;
          }
          dvadj[Fidx[fc - 1] + j - 1] = fc;
          free( f1 );
        }
      }

      curre = nexte;

      if ( Sixidx[curre - 1] == fin[k - 1] )
        nexte = Sixidx[curre - 1 + 2 * M_plan];
      else
        nexte = Sixidx[curre - 1 + 4 * M_plan];
    }

    Fidx[Nd] = Fidx[Nd - 1] + ( cnt );

    for ( i = 1; i <= cnt; i++ )
      dvadj[Fidx[Nd - 1] + i - 1] = newv[i - 1];

    int totaldone = 1;
    bool* done = (bool*)malloc( Nd * sizeof( bool ) );
    int jj;
    for ( jj = 0; jj < Nd; jj++ )
      done[jj] = false;

    done[Nd - 1] = true;
    int* looknext = (int*)malloc( ( Nd + 1 ) * sizeof( int ) );
    for ( jj = 0; jj < Nd + 1; jj++ )
      looknext[jj] = 0;
    looknext[0] = Nd;
    int looknextcount = 1;
    int index = 0;
    int count = 0;

    while ( looknext[index] != 0 )
    {
      int v = looknext[index];

      for ( i = Fidx[v - 1]; i < Fidx[v]; i++ )
      {
        int check = ( dvadj[i] );
        if ( done[check - 1] == false )
        {
          looknext[looknextcount] = check;
          looknextcount = looknextcount + 1;
          dist[k - 1 + ( check - 1 ) * d] = dist[k - 1 + ( v - 1 ) * d] + 1;
          prev[k - 1 + ( check - 1 ) * d] = v;
          done[check - 1] = true;
          totaldone = totaldone + 1;
        }
      }
      index = index + 1;
    }

    dist[k - 1 + ( Nd - 1 ) * d] = Nd;
    free( dvadj );
    free( Fidx );
    free( newv );
    free( done );
    free( looknext );
  }
  for ( i = 1; i <= d; i++ )
    dist[i - 1 + ( Nd - 1 ) * d] = Mdm;

  int* distsum = (int*)malloc( Nd * sizeof( int ) );
  for ( i = 1; i <= Nd; i++ )
    distsum[i - 1] = 0;

  for ( i = 1; i <= d; i++ )
    for ( j = 1; j <= Nd; j++ )
      distsum[j - 1] = distsum[j - 1] + dist[i - 1 + ( j - 1 ) * d];

  newcross = Mdm;
  while ( newface == 0 )
  {
    newcross = 2 * newcross;
    for ( i = 1; i <= Nd; i++ )
      if ( distsum[i - 1] < newcross )
      {
        newcross = distsum[i - 1];
        newface = i;
      }
  }

  improvementfound = 1;

  bool* checked = (bool*)malloc( Nd * sizeof( bool ) );
  for ( i = 0; i < Nd; i++ )
    checked[i] = false;
  maxpdist = 0;
  for ( i = 1; i <= d; i++ )
  {
    pdist[i - 1] = dist[i - 1 + ( newface - 1 ) * d];
    if ( pdist[i - 1] > maxpdist )
      maxpdist = pdist[i - 1];
  }

  pidx[0] = 0;
  pathsidx[0] = 0;
  for ( i = 1; i <= d; i++ )
  {
    pidx[i] = pidx[i - 1] + pdist[i - 1];
    pathsidx[i] = pathsidx[i - 1] + pdist[i - 1] + 1;
  }
  *paths = (int*)malloc( pathsidx[d] * sizeof( int ) );
  for ( i = 0; i < pathsidx[d]; i++ )
    ( *paths )[i] = 0;
  if ( maxpdist > 0 )
  {
    *pathedge = (int*)malloc( pidx[d] * sizeof( int ) );
    for ( i = 0; i < pidx[d]; i++ )
      ( *pathedge )[i] = 0;
  }
  else
  {
    *pathedge = (int*)malloc( d * sizeof( int ) );
    for ( i = 0; i < d; i++ )
      ( *pathedge )[i] = 0;
  }
  for ( k = 1; k <= d; k++ )
  {
    int nn = newface;
    ( *paths )[pathsidx[k - 1]] = nn;

    for ( j = 1; j <= pdist[k - 1] - 1; j++ )
    {
      nn = prev[k - 1 + ( nn - 1 ) * d];
      ( *paths )[pathsidx[k - 1] + j] = nn;
      if ( checked[nn - 1] )
        treecheck = 1;
      checked[nn - 1] = true;
    }

    if ( pdist[k - 1] > 0 )
    {
      int r;
      for ( r = 1; r <= maxdegfin; r++ )
      {
        if ( delv[k - 1 + ( r - 1 ) * d] == 0 )
          break;

        bool breakloop = false;

        for ( i = Fidxm[delv[k - 1 + ( r - 1 ) * d] - 1] + 1; i <= Fidxm[delv[k - 1 + ( r - 1 ) * d]]; i++ )
          if ( dvadjm[i - 1] == ( *paths )[pathsidx[k - 1] + ( pdist[k - 1] - 1 )] )
          {
            ( *paths )[pathsidx[k - 1] + ( pdist[k - 1] )] = delv[k - 1 + ( r - 1 ) * d];
            if ( checked[delv[k - 1 + ( r - 1 ) * d] - 1] )
              treecheck = 1;
            checked[delv[k - 1 + ( r - 1 ) * d] - 1] = true;
            ( *pathedge )[pidx[k - 1] + ( pdist[k - 1] - 1 )] = deadjm[i - 1];
            pathedgeonce[deadjm[i - 1] - 1] = k;
            pathedgeonce[deadjm[i - 1] - 1 + maxdeadjm] = pdist[k - 1];
            breakloop = true;
            break;
          }
        if ( breakloop )
          break;
      }
    }
  }

  if ( maxpdist > 0 )
  {
    for ( k = 1; k <= d; k++ )
    {
      for ( j = 1; j <= pdist[k - 1] - 1; j++ )
      {
        for ( i = Fidxm[( *paths )[pathsidx[k - 1] + ( j - 1 )] - 1] + 1; i <= Fidxm[( *paths )[pathsidx[k - 1] + ( j - 1 )]]; i++ )
        {
          if ( dvadjm[i - 1] == ( *paths )[pathsidx[k - 1] + j] )
          {
            ( *pathedge )[pidx[k - 1] + ( j - 1 )] = deadjm[i - 1];
            pathedgeonce[deadjm[i - 1] - 1] = k;
            pathedgeonce[deadjm[i - 1] - 1 + maxdeadjm] = j;
            break;
          }
        }
      }
    }
  }

  free( checked );
  free( dist );
  free( prev );
  free( distsum );

  int maxmaxpaths = 0;
  for ( i = 1; i <= pathsidx[d]; i++ )
    if ( ( *paths )[i - 1] > maxmaxpaths )
      maxmaxpaths = ( *paths )[i - 1];
  outputs[0] = newcross;
  outputs[1] = treecheck;
  outputs[2] = newface;
  outputs[3] = maxpdist;
  outputs[4] = improvementfound;
  outputs[5] = used_bigface;
  outputs[6] = maxmaxpaths;
}

void Planarise_Sixidx_PE( int* Sixidx, int* Clabel, int* Cidx, int* Dlabel, int N, int M, int* Sixidx_plan, int* planvlabels, int* elist, int* eidx, int* origedgelabel, int* origedgeorder, int* incidente_plan, int* newdim, int* not_added_edges )
{

  /******************************************************************************/
  /*                                                                            */
  /*  Planarise_Sixidx_PE                                                       */
  /*                                                                            */
  /*  A stripped down version of Planarise_Sixidx used in the planar embedding. */
  /*                                                                            */
  /******************************************************************************/

  int i;
  int j;
  int k;

  int* otherplanvlabel = (int*)malloc( Cidx[M] * sizeof( int ) );

  for ( i = 0; i < Cidx[M]; i++ )
  {
    planvlabels[i] = 0;
    otherplanvlabel[i] = 0;
  }
  for ( i = 0; i < M + Cidx[M]; i++ )
    elist[i] = 0;

  for ( i = 0; i < M + 1; i++ )
    eidx[i] = i + Cidx[i];

  for ( i = 0; i < M; i++ )
    elist[eidx[i]] = i + 1;

  int oldN = N;
  int oldM = M;

  for ( i = 0; i < M; i++ )
  {
    origedgelabel[i] = i + 1;
    origedgeorder[i] = 1;
  }
  for ( i = 0; i < Cidx[M]; i++ )
  {
    origedgelabel[M + i] = 0;
    origedgeorder[M + i] = 0;
  }

  for ( i = 0; i < N; i++ )
    incidente_plan[i] = 0;

  int nss = M + Cidx[M];
  for ( i = 0; i < M; i++ )
    for ( j = 0; j < 6; j++ )
      Sixidx_plan[i + j * nss] = Sixidx[i + j * M];

  for ( i = 0; i < Cidx[M]; i++ )
    for ( j = 0; j < 6; j++ )
      Sixidx_plan[M + i + j * nss] = 0;

  for ( i = 1; i <= oldM; i++ )
  {
    int tt = -1;
    if ( Cidx[i] - Cidx[i - 1] != 0 )
      tt = Sixidx_plan[( i - 1 ) + nss];

    for ( j = 1; j <= Cidx[i] - Cidx[i - 1]; j++ )
    {
      M = M + 1;
      elist[eidx[i - 1] + j] = M;
      origedgelabel[M - 1] = i;
      origedgeorder[M - 1] = j + 1;

      if ( planvlabels[Cidx[i - 1] + j - 1] == 0 )
      {
        N = N + 1;
        planvlabels[Cidx[i - 1] + j - 1] = N;
        int ce = Clabel[Cidx[i - 1] + j - 1];
        for ( k = Cidx[ce - 1] + 1; k <= Cidx[ce]; k++ )
        {
          if ( Clabel[k - 1] == i )
          {
            planvlabels[k - 1] = N;
            otherplanvlabel[Cidx[i - 1] + j - 1] = k;
            otherplanvlabel[k - 1] = Cidx[i - 1] + j;
          }
        }
      }
    }

    if ( Cidx[i] - Cidx[i - 1] > 0 )
    {
      Sixidx_plan[M - 1] = tt;
      Sixidx_plan[M - 1 + nss] = planvlabels[Cidx[i] - 1];
      Sixidx_plan[M - 1 + 2 * nss] = Sixidx_plan[i - 1 + 4 * nss];
      Sixidx_plan[M - 1 + 3 * nss] = Sixidx_plan[i - 1 + 5 * nss];
      Sixidx_plan[M - 1 + 4 * nss] = 0;
      Sixidx_plan[M - 1 + 5 * nss] = 0;
      Sixidx_plan[i - 1 + nss] = planvlabels[Cidx[i - 1]];
    }

    for ( j = 1; j <= Cidx[i] - Cidx[i - 1] - 1; j++ )
    {
      if ( planvlabels[Cidx[i - 1] + j - 1] < planvlabels[Cidx[i - 1] + j] )
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1] = planvlabels[Cidx[i - 1] + j - 1];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + nss] = planvlabels[Cidx[i - 1] + j];
      }
      else
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1] = planvlabels[Cidx[i - 1] + j];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + nss] = planvlabels[Cidx[i - 1] + j - 1];
      }
    }
  }

  for ( i = 1; i <= M; i++ )
  {
    if ( ( i <= oldM && not_added_edges[i - 1] == 0 ) || i > oldM )
    {
      if ( Sixidx_plan[i - 1] <= oldN )
      {
        int ec = Sixidx_plan[i - 1 + 2 * nss];
        if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1] )
          Sixidx_plan[i - 1 + 2 * nss] = elist[eidx[ec] - 1];

        ec = Sixidx_plan[i - 1 + 3 * nss];
        if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1] )
          Sixidx_plan[i - 1 + 3 * nss] = elist[eidx[ec] - 1];
      }
      if ( Sixidx_plan[i - 1 + nss] <= oldN )
      {
        int ec = Sixidx_plan[i - 1 + 4 * nss];
        if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1 + nss] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1 + nss] )
          Sixidx_plan[i - 1 + 4 * nss] = elist[eidx[ec] - 1];

        ec = Sixidx_plan[i - 1 + 5 * nss];
        if ( Sixidx_plan[ec - 1] != Sixidx_plan[i - 1 + nss] && Sixidx_plan[ec - 1 + nss] != Sixidx_plan[i - 1 + nss] )
          Sixidx_plan[i - 1 + 5 * nss] = elist[eidx[ec] - 1];
      }
    }
  }

  for ( i = 1; i <= oldM; i++ )
  {
    for ( j = 1; j <= Cidx[i] - Cidx[i - 1]; j++ )
    {
      int othere = Clabel[Cidx[i - 1] + j - 1];
      int f1 = otherplanvlabel[Cidx[i - 1] + j - 1] - Cidx[othere - 1];
      int a1 = 1;
      int a2 = 0;

      if ( Dlabel[Cidx[i - 1] + j - 1] != 1 )
      {
        a1 = 0;
        a2 = 1;
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1] )
      {
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 2 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 3 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
      }
      else
      {
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 4 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
        Sixidx_plan[elist[eidx[i - 1] + j - 1] - 1 + 5 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[i - 1] + j] - 1] )
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 2 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 3 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
      }
      else
      {
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 4 * nss] = elist[eidx[othere - 1] + f1 - 1 + a2];
        Sixidx_plan[elist[eidx[i - 1] + j] - 1 + 5 * nss] = elist[eidx[othere - 1] + f1 - 1 + a1];
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1] )
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 2 * nss] = elist[eidx[i - 1] + j - 1 + a2];
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 3 * nss] = elist[eidx[i - 1] + j - 1 + a1];
      }
      else
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 4 * nss] = elist[eidx[i - 1] + j - 1 + a2];
        Sixidx_plan[elist[eidx[othere - 1] + f1 - 1] - 1 + 5 * nss] = elist[eidx[i - 1] + j - 1 + a1];
      }

      if ( planvlabels[Cidx[i - 1] + j - 1] == Sixidx_plan[elist[eidx[othere - 1] + f1] - 1] )
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 2 * nss] = elist[eidx[i - 1] + j - 1 + a1];
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 3 * nss] = elist[eidx[i - 1] + j - 1 + a2];
      }
      else
      {
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 4 * nss] = elist[eidx[i - 1] + j - 1 + a1];
        Sixidx_plan[elist[eidx[othere - 1] + f1] - 1 + 5 * nss] = elist[eidx[i - 1] + j - 1 + a2];
      }
    }
  }

  for ( i = 1; i <= M; i++ )
  {
    if ( ( i <= oldM && not_added_edges[i - 1] == 0 ) || ( i > oldM ) )
    {
      if ( Sixidx_plan[i - 1] <= oldN && incidente_plan[Sixidx_plan[i - 1] - 1] == 0 )
        incidente_plan[Sixidx_plan[i - 1] - 1] = i;
      if ( Sixidx_plan[i - 1 + nss] <= oldN && incidente_plan[Sixidx_plan[i - 1 + nss] - 1] == 0 )
        incidente_plan[Sixidx_plan[i - 1 + nss] - 1] = i;
    }
  }

  newdim[0] = N;
  newdim[1] = M;

  free( otherplanvlabel );
}

void Embed_Planar_Subgraph( int** Sixidx, int** incidente, int** Clabel, int** Dlabel, int** Cidx, int N, int M, int** deg, int* vadj, int* vidx, int* eadj, int seed, int maxdeg, int* outputs, int verbose )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Embed_Planar_Subgraph                                                  */
  /*                                                                         */
  /*  Produces a valid embedding of the graph by the following method:       */
  /*                                                                         */
  /*  First, a small cycle is identified with no inner edges. The cycle is a */
  /*  vertex-induced subgraph. If none exists then the graph is acyclic and  */
  /*  the crossing number is 0.                                              */
  /*                                                                         */
  /*  Then, one at a time, vertices not in the subgraph are reintroduced.    */
  /*  A vertex will only be reintroduced if it has at least one neighbour    */
  /*  in the subgraph. The vertex is placed in the optimal face using the    */
  /*  same approach as in Quick_Cross_Main_Loop, except a few features are   */
  /*  not required. For example, Undo_Sub is not needed because only new     */
  /*  vertices are added each iteration, so subdivisions remain necessary.   */
  /*                                                                         */
  /*  Embed_Planar_Subgraph calls functions from the main code, as well as   */
  /*  stripped down versions of Planarise_Sixidx, Get_Distances, New_Clabel  */
  /*  and Subd_Try.                                                          */
  /*                                                                         */
  /***************************************************************************/

  int* dist = (int*)malloc( N * sizeof( int ) );
  int* done = (int*)malloc( N * sizeof( int ) );
  int i;
  int* p = (int*)malloc( N * sizeof( int ) );
  int* looknext = (int*)malloc( ( N + 1 ) * sizeof( int ) );
  looknext[N] = 0;
  for ( i = 0; i < N; i++ )
  {
    dist[i] = 0;
    done[i] = 0;
    p[i] = 0;
    looknext[i] = 0;
  }

  int totaldone = 1;
  int looknextcount = 1;
  int index = 1;
  int v = ( rand() % N ) + 1;
  int t = 0;
  int isave;

  done[v - 1] = 1;
  looknext[0] = v;
  while ( looknext[index - 1] != 0 )
  {
    v = looknext[index - 1];
    for ( i = vidx[v - 1]; i < vidx[v]; i++ )
    {
      isave = i;
      if ( done[vadj[i] - 1] == 0 )
      {
        looknextcount = looknextcount + 1;
        looknext[looknextcount - 1] = vadj[i];
        dist[vadj[i] - 1] = dist[v - 1] + 1;
        p[vadj[i] - 1] = v;
        done[vadj[i] - 1] = 1;
        totaldone = totaldone + 1;
      }
      else
      {
        if ( p[v - 1] != vadj[i] )
        {
          t = 1;
          break;
        }
      }
    }
    if ( t > 0 )
    {
      break;
    }
    index = index + 1;
  }
  if ( t == 0 )
  { // then acyclic
    return;
  }
  int otherv = vadj[isave];
  int dv1 = dist[v - 1];
  int dv2 = dist[otherv - 1];
  int* not_added_edges = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
  {
    not_added_edges[i] = 1;
  }
  int* checked = (int*)malloc( N * sizeof( int ) );
  int j;
  int ef;
  for ( i = 0; i < N; i++ )
  {
    checked[i] = 0;
  }
  int eb1;
  for ( i = vidx[v - 1]; i < vidx[v]; i++ )
  {
    if ( vadj[i] == otherv )
    {
      eb1 = eadj[i];
      break;
    }
  }
  checked[v - 1] = 1;
  checked[otherv - 1] = 1;
  int eb2 = eb1;
  int swt;
  if ( dv1 > dv2 )
  {
    swt = 0;
  }
  else
  {
    swt = 1;
  }

  not_added_edges[eb1 - 1] = 0;
  int vb;
  int ef1;
  int ef2;

  for ( i = 0; i < dv1 + dv2 + 1; i++ )
  {
    if ( swt == 0 )
    {
      vb = v;
      v = p[v - 1];
      swt = 1;
      for ( j = vidx[v - 1]; j < vidx[v]; j++ )
      {
        if ( vadj[j] == vb )
        {
          ef1 = eadj[j];
          break;
        }
      }
      not_added_edges[ef1 - 1] = 0;
      if ( ( *Sixidx )[eb1 - 1] == vb )
      {
        ( *Sixidx )[2 * M + eb1 - 1] = ef1;
        ( *Sixidx )[3 * M + eb1 - 1] = ef1;
      }
      else
      {
        ( *Sixidx )[4 * M + eb1 - 1] = ef1;
        ( *Sixidx )[5 * M + eb1 - 1] = ef1;
      }
      if ( ( *Sixidx )[ef1 - 1] == vb )
      {
        ( *Sixidx )[2 * M + ef1 - 1] = eb1;
        ( *Sixidx )[3 * M + ef1 - 1] = eb1;
      }
      else
      {
        ( *Sixidx )[4 * M + ef1 - 1] = eb1;
        ( *Sixidx )[5 * M + ef1 - 1] = eb1;
      }
      eb1 = ef1;
      if ( checked[v - 1] == 0 )
      {
        checked[v - 1] = 1;
      }
      else
      {
        break;
      }
    }
    else
    {
      vb = otherv;
      otherv = p[otherv - 1];
      swt = 0;
      for ( j = vidx[otherv - 1]; j < vidx[otherv]; j++ )
      {
        if ( vadj[j] == vb )
        {
          ef2 = eadj[j];
          break;
        }
      }
      not_added_edges[ef2 - 1] = 0;
      if ( ( *Sixidx )[eb2 - 1] == vb )
      {
        ( *Sixidx )[2 * M + eb2 - 1] = ef2;
        ( *Sixidx )[3 * M + eb2 - 1] = ef2;
      }
      else
      {
        ( *Sixidx )[4 * M + eb2 - 1] = ef2;
        ( *Sixidx )[5 * M + eb2 - 1] = ef2;
      }
      if ( ( *Sixidx )[ef2 - 1] == vb )
      {
        ( *Sixidx )[2 * M + ef2 - 1] = eb2;
        ( *Sixidx )[3 * M + ef2 - 1] = eb2;
      }
      else
      {
        ( *Sixidx )[4 * M + ef2 - 1] = eb2;
        ( *Sixidx )[5 * M + ef2 - 1] = eb2;
      }
      eb2 = ef2;
      if ( checked[otherv - 1] == 0 )
      {
        checked[otherv - 1] = 1;
      }
      else
      {
        break;
      }
    }
  }
  if ( ( *Sixidx )[ef1 - 1] == v )
  {
    ( *Sixidx )[2 * M + ef1 - 1] = ef2;
    ( *Sixidx )[3 * M + ef1 - 1] = ef2;
  }
  else
  {
    ( *Sixidx )[4 * M + ef1 - 1] = ef2;
    ( *Sixidx )[5 * M + ef1 - 1] = ef2;
  }
  if ( ( *Sixidx )[ef2 - 1] == v )
  {
    ( *Sixidx )[2 * M + ef2 - 1] = ef1;
    ( *Sixidx )[3 * M + ef2 - 1] = ef1;
  }
  else
  {
    ( *Sixidx )[4 * M + ef2 - 1] = ef1;
    ( *Sixidx )[5 * M + ef2 - 1] = ef1;
  }

  if ( verbose > 0 )
    fprintf( stdout, "\nStarting planar embedding with %d vertices and 0 crossings.\n", dv1 + dv2 + 1 );

  int current_crossings = 0;
  int vm = 0;
  int* verts; // TO DO - see how many subdivisions occur and rescale this appropriately.
  verts = (int*)malloc( N * sizeof( int ) );
  int d;
  int* fin = (int*)malloc( maxdeg * sizeof( int ) );
  int* fi = (int*)malloc( maxdeg * sizeof( int ) );
  int* fitemp = (int*)malloc( M * sizeof( int ) );
  int ficnt = 0;
  int outputssubd[3];
  int gdoutputs[7];
  int maxdeadjm;
  int maxdegfin;
  int newcross;
  int treecheck;
  int newface;
  int maxpdist;
  int maxmaxpaths;
  int Moriginal = M;
  int origN = N;
  int maxcrs = N;

  for ( i = 0; i < maxcrs; i++ )
  {
    ( *Clabel )[i] = 0;
    ( *Dlabel )[i] = 0;
  }

  while ( true )
  {
    vm = 0;

    if ( seed > 0 )
    {
      randperm( &verts, origN );
    }
    else
    {
      for ( i = 0; i < origN; i++ )
      {
        verts[i] = i + 1;
      }
    }
    for ( i = 0; i < origN; i++ )
    { // if we subdivide, new verts and edges are marked as inserted.
      if ( checked[verts[i] - 1] == 1 )
      {
        for ( j = vidx[verts[i] - 1]; j < vidx[verts[i]]; j++ )
        {
          if ( not_added_edges[eadj[j] - 1] == 1 )
          {
            if ( ( *Sixidx )[eadj[j] - 1] == verts[i] )
            {
              vm = ( *Sixidx )[M + eadj[j] - 1];
            }
            else
            {
              vm = ( *Sixidx )[eadj[j] - 1];
            }
            checked[vm - 1] = 1;
            break;
          }
        }
        if ( vm > 0 )
        {
          break;
        }
      }
    }
    if ( vm == 0 )
    {
      break;
    }

    int* Sixidx_plan = (int*)malloc( ( 6 * M + 6 * ( *Cidx )[M] ) * sizeof( int ) );
    int* planvlabels = (int*)malloc( ( ( *Cidx )[M] ) * sizeof( int ) );
    int* elist = (int*)malloc( ( M + ( *Cidx )[M] ) * sizeof( int ) );
    int* eidx = (int*)malloc( ( M + 1 ) * sizeof( int ) );
    int* origedgelabel = (int*)malloc( ( M + ( *Cidx )[M] ) * sizeof( int ) );
    int* origedgeorder = (int*)malloc( ( M + ( *Cidx )[M] ) * sizeof( int ) );
    int* incidente_plan = (int*)malloc( ( N ) * sizeof( int ) );
    int newdim[2];
    Planarise_Sixidx_PE( *Sixidx, *Clabel, *Cidx, *Dlabel, N, M, Sixidx_plan, planvlabels, elist, eidx, origedgelabel, origedgeorder, incidente_plan, newdim, not_added_edges );
    int N_plan = newdim[0];
    int M_plan = newdim[1];
    // faces uses original faces code, with fi as find(not_added edges)
    ficnt = -1;
    for ( i = 0; i < M; i++ )
    {
      if ( not_added_edges[i] == 1 )
      {
        ficnt = ficnt + 1;
        fitemp[ficnt] = i + 1;
      }
    }
    // faces stuff
    int* F = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* Fidx = (int*)malloc( ( ( 2 * M_plan + 5 ) / 3 ) * sizeof( int ) );
    int* Fdir = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* dvadj = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* ddeg = (int*)malloc( ( 2 + M_plan - N_plan ) * sizeof( int ) );
    int* dualel = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int* deadj = (int*)malloc( ( 2 * M_plan ) * sizeof( int ) );
    int Nd = Faces( Sixidx_plan, M_plan, N_plan, fitemp, ficnt + 1, F, Fidx, Fdir, dvadj, ddeg, dualel, deadj );
    d = 0;
    maxdegfin = 0;
    for ( i = vidx[vm - 1]; i < vidx[vm]; i++ )
    {
      if ( checked[vadj[i] - 1] == 1 )
      {
        d = d + 1;
        fin[d - 1] = vadj[i];
        fi[d - 1] = eadj[i];
        if ( ( *deg )[vadj[i] - 1] > maxdegfin )
        {
          maxdegfin = ( *deg )[vadj[i] - 1];
        }
      }
    }
    maxdeadjm = M_plan;

    int* pathedgeonce = (int*)malloc( ( 2 * maxdeadjm ) * sizeof( int ) );
    int* pdist = (int*)malloc( d * sizeof( int ) );
    int* pidx = (int*)malloc( ( d + 1 ) * sizeof( int ) );
    int* pathsidx = (int*)malloc( ( d + 1 ) * sizeof( int ) );
    int* delv = (int*)malloc( d * maxdegfin * sizeof( int ) );
    int gdoutputs[7];
    int* bad_edges = (int*)malloc( M * sizeof( int ) );
    for ( i = 0; i < M; i++ )
    {
      bad_edges[i] = 0;
    }
    for ( i = 0; i < 2 * maxdeadjm; i++ )
    {
      pathedgeonce[i] = 0;
    }
    for ( i = 0; i < d * maxdegfin; i++ )
    {
      delv[i] = 0;
    }

    // getdist
    int* pathedge;
    int* paths;
    Get_Distances_PE( dvadj, Fidx, Nd, M_plan, M_plan, Sixidx_plan, dualel, 0, d, fin, incidente_plan, ( *deg ), current_crossings, deadj, maxdegfin, maxdeadjm, &pathedge, pidx, pathedgeonce, &paths, pathsidx, pdist, delv, gdoutputs );
    newcross = gdoutputs[0];
    treecheck = gdoutputs[1];
    newface = gdoutputs[2];
    maxpdist = gdoutputs[3];
    maxmaxpaths = gdoutputs[6];
    current_crossings = current_crossings + newcross;
    if ( 2 * ( newcross + current_crossings ) > maxcrs )
    {
      // then reallocate Clab and Dlab vectors
      maxcrs = 4 * ( newcross + current_crossings );
      *Clabel = (int*)realloc( *Clabel, maxcrs * sizeof( int ) );
      *Dlabel = (int*)realloc( *Dlabel, maxcrs * sizeof( int ) );
    }
    New_Clabel_PE( M, M_plan, Sixidx_plan, ( *Sixidx ), ( *Clabel ), ( *Dlabel ), ( *Cidx ), elist, eidx, pathedge, paths, pidx, pathsidx, pdist, maxpdist, d, fin, fi, F, Fdir, Fidx, vm, origedgelabel, origedgeorder, newface, pathedgeonce, treecheck, maxdeadjm );

    for ( i = 0; i < d; i++ )
    {
      not_added_edges[fi[i] - 1] = 0;
    }
    Subd_Try_PE( current_crossings, Sixidx, ( *Clabel ), ( *Dlabel ), Cidx, pathedge, N, M, origedgelabel, incidente, deg, pdist, fi, fin, vm, d, maxpdist, pidx, outputssubd, &not_added_edges, &checked );
    N = outputssubd[0];
    M = outputssubd[1];
    current_crossings = outputssubd[2];

    free( Sixidx_plan );
    free( planvlabels );
    free( elist );
    free( eidx );
    free( origedgelabel );
    free( origedgeorder );
    free( incidente_plan );
    free( F );
    free( Fidx );
    free( Fdir );
    free( dvadj );
    free( ddeg );
    free( dualel );
    free( deadj );
    free( pathedgeonce );
    free( pdist );
    free( pidx );
    free( pathsidx );
    free( pathedge );
    free( paths );
    free( delv );
    free( bad_edges );

    if ( verbose > 0 )
      fprintf( stdout, "Planar Embedding: Added vertex: %d, current crossings: %d\n", vm, current_crossings );
  }

  free( dist );
  free( done );
  free( p );
  free( looknext );
  free( not_added_edges );
  free( checked );
  free( verts );
  free( fin );
  free( fi );
  free( fitemp );

  outputs[0] = N;
  outputs[1] = M;
  outputs[2] = current_crossings;
}

void Embed_Circ( int** Clabel, int** Dlabel, int* Cidx, int* Sixidx, int N, int M, int maxdeg, int* vadj, int* eadj, int* vidx, int seed )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Embed_Circ                                                             */
  /*                                                                         */
  /*  Produces a valid embedding of the graph by the following method:       */
  /*                                                                         */
  /*  Each vertex is placed on a circle. This ensures that all edges are     */
  /*  inside the circle, and it is impossible for any two intersecting edges */
  /*  to be parallel. To avoid vertical edges (to make calculations simpler) */
  /*  a dummy vertex is added is the number of vertices is even.             */
  /*                                                                         */
  /***************************************************************************/

  double pi = 3.14159265358979323846;
  int* smap;
  int* invsmap;
  int i;
  int N0;
  smap = (int*)malloc( N * sizeof( int ) );
  invsmap = (int*)malloc( N * sizeof( int ) );
  if ( seed > 0 )
  {
    randperm( &smap, N );
    for ( i = 0; i < N; i++ )
    {
      invsmap[smap[i] - 1] = i + 1;
    }
  }
  else
  {
    for ( i = 0; i < N; i++ )
    {
      smap[i] = i + 1;
      invsmap[i] = i + 1;
    }
  }

  double* cv2 = (double*)malloc( 2 * N * sizeof( double ) );
  ;
  for ( i = 0; i < 2 * N; i++ )
  {
    cv2[i] = 0;
  }
  if ( N % 2 > 0 )
  {
    for ( i = 0; i < N; i++ )
    {
      cv2[i] = sin( ( 2 * pi * ( i + 1 ) ) / N );
      cv2[N + i] = cos( ( 2 * pi * ( i + 1 ) ) / N );
    }
    N0 = N;
  }
  else
  {
    N0 = N + 1;
    for ( i = 0; i < N; i++ )
    {
      cv2[i] = sin( ( 2 * pi * ( i + 1 ) ) / ( N0 ) );
      cv2[N + i] = cos( ( 2 * pi * ( i + 1 ) ) / ( N0 ) );
    }
  }

  double* cv = (double*)malloc( 2 * N * sizeof( double ) );
  for ( i = 0; i < N; i++ )
  {
    cv[i] = cv2[smap[i] - 1];
    cv[N + i] = cv2[N + smap[i] - 1];
  }

  int fpos;
  int startpos;
  int f1;
  int k;
  int fs;
  int curre;
  int starte;
  int j;
  int fsum;
  int currpos;

  for ( i = 0; i < N; i++ )
  {
    startpos = smap[i];
    fpos = startpos;
    f1 = -1;
    while ( f1 == -1 )
    {
      if ( fpos == N )
      {
        fpos = 1;
      }
      else
      {
        fpos = fpos + 1;
      }
      for ( k = 0; k < N; k++ )
      {
        if ( smap[k] == fpos )
        {
          fs = k + 1;
          break;
        }
      }

      for ( k = vidx[i]; k < vidx[i + 1]; k++ )
      {
        if ( vadj[k] == fs )
        {
          f1 = k;
          break;
        }
      }
    }

    currpos = startpos;
    curre = eadj[f1];
    starte = curre;
    for ( k = 0; k < N - 1; k++ )
    {
      if ( currpos < N )
      {
        currpos = currpos + 1;
      }
      else
      {
        currpos = 1;
      }
      if ( currpos == startpos )
      {
        break;
      }
      for ( j = 0; j < N; j++ )
      {
        if ( smap[j] == currpos )
        {
          fs = j + 1;
          break;
        }
      }
      fsum = 0;
      for ( j = vidx[i]; j < vidx[i + 1]; j++ )
      {
        if ( vadj[j] == fs )
        {
          fsum = 1;
          f1 = j;
          break;
        }
      }
      if ( fsum == 1 )
      {
        if ( Sixidx[eadj[f1] - 1] == i + 1 )
        {
          Sixidx[3 * M + eadj[f1] - 1] = curre;
        }
        else
        {
          Sixidx[5 * M + eadj[f1] - 1] = curre;
        }
        if ( Sixidx[curre - 1] == i + 1 )
        {
          Sixidx[2 * M + curre - 1] = eadj[f1];
        }
        else
        {
          Sixidx[4 * M + curre - 1] = eadj[f1];
        }
        curre = eadj[f1];
      }
    }
    if ( Sixidx[curre - 1] == i + 1 )
    {
      Sixidx[2 * M + curre - 1] = starte;
    }
    else
    {
      Sixidx[4 * M + curre - 1] = starte;
    }
    if ( Sixidx[starte - 1] == i + 1 )
    {
      Sixidx[3 * M + starte - 1] = curre;
    }
    else
    {
      Sixidx[5 * M + starte - 1] = curre;
    }
  }
  int maxcrs = N;
  int* first;
  int* next;
  int* crsdir;
  int* crslabel;
  double* crscoord;
  int mcnt = 0;
  first = (int*)calloc( M, sizeof( int ) );
  next = (int*)calloc( maxcrs, sizeof( int ) );
  crscoord = (double*)calloc( maxcrs, sizeof( double ) );
  crslabel = (int*)calloc( maxcrs, sizeof( int ) );
  crsdir = (int*)calloc( maxcrs, sizeof( int ) );
  int swap;
  int swap2;
  int* H;
  int* L;
  H = (int*)malloc( N * sizeof( int ) );
  L = (int*)malloc( N * sizeof( int ) );
  int Hcnt;
  int Lcnt;
  int prev;
  int prev2;
  int e1;
  int v1;
  double cos1;
  double cos2;
  double xcross;
  int extreme = 0;
  int a1;
  int a2;
  int b1;
  int b2;
  int r;
  int jj;
  int kk;

  for ( i = 1; i <= M; i++ )
  {
    a1 = smap[Sixidx[i - 1] - 1];
    a2 = smap[Sixidx[i - 1 + M] - 1];
    swap2 = 0;
    Hcnt = 0;
    Lcnt = 0;

    if ( a2 > a1 )
    {
      swap = 0;
      for ( j = a1 - 1; j >= 1; j-- )
      {
        Hcnt++;
        H[Hcnt - 1] = j;
      }
      for ( j = N; j >= a2 + 1; j-- )
      {
        Hcnt++;
        H[Hcnt - 1] = j;
      }
      for ( j = a1 + 1; j <= a2 - 1; j++ )
      {
        Lcnt++;
        L[Lcnt - 1] = j;
      }
      if ( cv[Sixidx[i - 1] - 1] > cv[Sixidx[i - 1 + M] - 1] )
      {
        swap2 = 1;
      }
    }
    else
    {
      swap = 1;
      for ( j = a1 + 1; j <= N; j++ )
      {
        Hcnt++;
        H[Hcnt - 1] = j;
      }
      for ( j = 1; j <= a2 - 1; j++ )
      {
        Hcnt++;
        H[Hcnt - 1] = j;
      }
      for ( j = a1 - 1; j >= a2 + 1; j-- )
      {
        Lcnt++;
        L[Lcnt - 1] = j;
      }
      if ( cv[Sixidx[i - 1 + M] - 1] < cv[Sixidx[i - 1] - 1] )
      {
        swap2 = 1;
      }
    }

    for ( j = 1; j <= Lcnt; j++ )
    {
      jj = L[j - 1];
      prev = 0;
      prev2 = 0;
      e1 = 0;
      for ( k = 1; k <= Hcnt; k++ )
      {
        kk = H[k - 1];
        for ( r = vidx[invsmap[jj - 1] - 1]; r < vidx[invsmap[jj - 1]]; r++ )
        {
          if ( vadj[r] == invsmap[kk - 1] )
          {
            e1 = eadj[r];
            v1 = kk;
            break;
          }
        }
        if ( e1 > 0 )
        {
          break;
        }
      }
      if ( e1 > 0 )
      {
        cos1 = cos( pi / N0 * ( a2 - a1 ) );
        cos2 = cos( pi / N0 * ( a2 + a1 ) );

        for ( k = 1; k <= vidx[invsmap[jj - 1]] - vidx[invsmap[jj - 1] - 1]; k++ )
        {
          b1 = smap[Sixidx[e1 - 1] - 1];
          b2 = smap[Sixidx[e1 - 1 + M] - 1];

          xcross = ( cos1 * cos( pi / N0 * ( b1 + b2 ) ) - cos( pi / N0 * ( b2 - b1 ) ) * cos2 ) / sin( pi / N0 * ( a1 + a2 - b1 - b2 ) );
          mcnt++;
          if ( mcnt > maxcrs )
          {

            next = (int*)realloc( next, 4 * maxcrs * sizeof( int ) );
            crscoord = (double*)realloc( crscoord, 4 * maxcrs * sizeof( double ) );
            crslabel = (int*)realloc( crslabel, 4 * maxcrs * sizeof( int ) );
            crsdir = (int*)realloc( crsdir, 4 * maxcrs * sizeof( int ) );
            for ( r = maxcrs; r < 4 * maxcrs; r++ )
            {
              next[r] = 0;
            }
            maxcrs = 4 * maxcrs;
          }
          crscoord[mcnt - 1] = xcross;
          crslabel[mcnt - 1] = e1;

          if ( invsmap[v1 - 1] < invsmap[jj - 1] )
          {
            if ( swap > 0 )
            {
              crsdir[mcnt - 1] = -1;
            }
            else
            {
              crsdir[mcnt - 1] = 1;
            }
          }
          else
          {
            if ( swap > 0 )
            {
              crsdir[mcnt - 1] = 1;
            }
            else
            {
              crsdir[mcnt - 1] = -1;
            }
          }
          if ( first[i - 1] == 0 )
          {
            first[i - 1] = mcnt;
          }
          if ( prev == 0 )
          {
            prev = first[i - 1];
          }

          for ( r = 1; r <= M; r++ )
          {
            if ( swap2 == 0 )
            {
              if ( xcross > crscoord[prev - 1] )
              {
                extreme = 1;
              }
            }
            else
            {
              if ( xcross < crscoord[prev - 1] )
              {
                extreme = 1;
              }
            }

            if ( extreme > 0 )
            {
              extreme = 0;
              if ( next[prev - 1] != 0 )
              {
                prev2 = prev;
                prev = next[prev - 1];
              }
              else
              {
                next[prev - 1] = mcnt;
                break;
              }
              continue;
            }
            else
            {
              if ( prev != first[i - 1] )
              {
                next[mcnt - 1] = prev;
                next[prev2 - 1] = mcnt;
                prev = mcnt;
              }
              else
              {
                if ( first[i - 1] != mcnt )
                {
                  next[mcnt - 1] = first[i - 1];
                }
                first[i - 1] = mcnt;
                prev = mcnt;
              }
              break;
            }
          }

          if ( swap == 0 )
          {
            if ( Sixidx[e1 - 1] == invsmap[jj - 1] )
            {
              e1 = Sixidx[e1 - 1 + 3 * M];
            }
            else
            {
              e1 = Sixidx[e1 - 1 + 5 * M];
            }

            if ( Sixidx[e1 - 1] == invsmap[jj - 1] )
            {
              if ( smap[Sixidx[e1 - 1 + M] - 1] <= a2 && smap[Sixidx[e1 - 1 + M] - 1] >= a1 )
              {
                break;
              }
              v1 = smap[Sixidx[e1 - 1 + M] - 1];
            }
            else
            {
              if ( smap[Sixidx[e1 - 1] - 1] <= a2 && smap[Sixidx[e1 - 1] - 1] >= a1 )
              {
                break;
              }
              v1 = smap[Sixidx[e1 - 1] - 1];
            }
          }
          else
          {
            if ( Sixidx[e1 - 1] == invsmap[jj - 1] )
            {
              e1 = Sixidx[e1 - 1 + 2 * M];
            }
            else
            {
              e1 = Sixidx[e1 - 1 + 4 * M];
            }

            if ( Sixidx[e1 - 1] == invsmap[jj - 1] )
            {
              if ( smap[Sixidx[e1 - 1 + M] - 1] >= a2 && smap[Sixidx[e1 - 1 + M] - 1] <= a1 )
              {
                break;
              }
              v1 = smap[Sixidx[e1 - 1 + M] - 1];
            }
            else
            {
              if ( smap[Sixidx[e1 - 1] - 1] >= a2 && smap[Sixidx[e1 - 1] - 1] <= a1 )
              {
                break;
              }
              v1 = smap[Sixidx[e1 - 1] - 1];
            }
          }
        }
      }
    }
  }
  if ( mcnt > N )
  {
    *Clabel = (int*)realloc( *Clabel, mcnt * sizeof( int ) );
    *Dlabel = (int*)realloc( *Dlabel, mcnt * sizeof( int ) );
  }
  mcnt = -1;
  for ( i = 1; i <= M; i++ )
  {
    curre = first[i - 1];
    if ( curre != 0 )
    {
      mcnt++;
      ( *Clabel )[mcnt] = crslabel[curre - 1];
      ( *Dlabel )[mcnt] = crsdir[curre - 1];
      for ( j = 1; j <= M; j++ )
      {
        curre = next[curre - 1];
        if ( curre != 0 )
        {
          mcnt++;
          ( *Clabel )[mcnt] = crslabel[curre - 1];
          ( *Dlabel )[mcnt] = crsdir[curre - 1];
        }
        else
        {
          break;
        }
      }
    }
    Cidx[i] = mcnt + 1;
  }
  free( smap );
  free( cv2 );
  free( cv );
  free( H );
  free( L );
  free( first );
  free( next );
  free( crscoord );
  free( crslabel );
  free( invsmap );
  free( crsdir );
}

void KKSL( long double* x, long double* y, int* Sixidx, int N, int* deg, int* vadj, int* vidx, int k, double l, int seed )
{

  /***************************************************************************/
  /*                                                                         */
  /*  KKSL                                                                   */
  /*                                                                         */
  /*  An implementation of the Kamada Kawai Spring Layout method.            */
  /*                                                                         */
  /***************************************************************************/

  long double tolerance = 0.001 / N;
  long double pi = acosl( -1.0L );
  int i;
  int j;
  // ell's here not 1's
  int* smap;
  // int *invsmap;
  smap = (int*)malloc( N * sizeof( int ) );
  // invsmap = (int *)malloc(N*sizeof(int));
  if ( seed > 0 )
  {
    randperm( &smap, N );
    // for(i=0; i<N; i++){
    // invsmap[smap[i]-1] = i+1;
    //}
  }
  else
  {
    for ( i = 0; i < N; i++ )
    {
      smap[i] = i + 1;
      // invsmap[i] = i+1;
    }
  }
  for ( i = 1; i <= N; i++ )
  {
    long double tt = 2 * pi * smap[i - 1] / N;

    x[i - 1] = l / 2.0 + l / 2.0 * cosl( tt );
    y[i - 1] = l / 2.0 + l / 2.0 * sinl( tt );
  }
  int* D;
  D = (int*)calloc( N * N, sizeof( int ) );
  int* stack;
  int* checked;
  stack = (int*)calloc( N, sizeof( int ) );
  checked = (int*)calloc( N, sizeof( int ) );
  int stackcurr = 0;
  int stackmax = 0;
  int currv;
  for ( i = 1; i <= N; i++ )
  {
    stack[0] = i;
    stackcurr = -1;
    stackmax = 0;
    for ( j = 1; j <= N; j++ )
    {
      checked[j - 1] = 0;
    }
    checked[i - 1] = 1;
    while ( stackcurr < stackmax )
    {
      stackcurr++;
      currv = stack[stackcurr];
      for ( j = vidx[currv - 1]; j < vidx[currv]; j++ )
      {
        if ( checked[vadj[j] - 1] == 0 )
        {
          D[N * ( i - 1 ) + vadj[j] - 1] = D[N * ( i - 1 ) + currv - 1] + 1;
          checked[vadj[j] - 1] = 1;
          stackmax++;
          stack[stackmax] = vadj[j];
        }
      }
    }
  }
  int mmD = 0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = i + 1; j < N; j++ )
    {
      if ( D[N * j + i] > mmD )
      {
        mmD = D[N * j + i];
      }
    }
  }
  long double* L;
  L = (long double*)malloc( N * N * sizeof( long double ) );
  long double* K;
  K = (long double*)malloc( N * N * sizeof( long double ) );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      L[N * ( j - 1 ) + i - 1] = ( (long double)l / (long double)mmD ) * (long double)D[N * ( j - 1 ) + i - 1];
      K[N * ( j - 1 ) + i - 1] = (long double)k / ( (long double)D[N * ( j - 1 ) + i - 1] * (long double)D[N * ( j - 1 ) + i - 1] );
    }
  }
  long double* T;
  T = (long double*)calloc( N * N, sizeof( long double ) );
  for ( i = 1; i <= N; i++ )
  {
    for ( j = i + 1; j <= N; j++ )
    {
      T[N * ( j - 1 ) + i - 1] = powl( ( x[j - 1] - x[i - 1] ), 2 ) + powl( ( y[j - 1] - y[i - 1] ), 2 );
      T[N * ( i - 1 ) + j - 1] = T[N * ( j - 1 ) + i - 1];
    }
  }
  long double* par_x;
  par_x = (long double*)malloc( N * sizeof( long double ) );
  long double* par_y;
  par_y = (long double*)malloc( N * sizeof( long double ) );
  for ( i = 1; i <= N; i++ )
  {
    par_x[i - 1] = 0;
    par_y[i - 1] = 0;
  }

  long double* Delta;
  Delta = (long double*)malloc( N * sizeof( long double ) );
  int m;
  for ( m = 1; m <= N; m++ )
  {
    for ( i = 1; i <= N; i++ )
    {
      if ( i == m )
      {
        continue;
      }
      par_x[m - 1] = par_x[m - 1] + K[( i - 1 ) * N + m - 1] * ( ( x[m - 1] - x[i - 1] ) - (long double)L[( i - 1 ) * N + m - 1] * ( x[m - 1] - x[i - 1] ) / sqrtl( powl( ( x[m - 1] - x[i - 1] ), 2 ) + powl( ( y[m - 1] - y[i - 1] ), 2 ) ) );
      par_y[m - 1] = par_y[m - 1] + K[( i - 1 ) * N + m - 1] * ( ( y[m - 1] - y[i - 1] ) - (long double)L[( i - 1 ) * N + m - 1] * ( y[m - 1] - y[i - 1] ) / sqrtl( powl( ( x[m - 1] - x[i - 1] ), 2 ) + powl( ( y[m - 1] - y[i - 1] ), 2 ) ) );
    }
    Delta[m - 1] = sqrtl( powl( par_x[m - 1], 2 ) + powl( par_y[m - 1], 2 ) );
  }

  long double a = -1000000;
  long double prev_max_delta = 0;

  for ( i = 0; i < N; i++ )
  {
    if ( Delta[i] > a )
    {
      a = Delta[i];
      m = i + 1;
    }
  }

  long double* contributionx_m;
  long double* contributiony_m;
  contributionx_m = (long double*)malloc( N * sizeof( long double ) );
  contributiony_m = (long double*)malloc( N * sizeof( long double ) );
  int first_time;
  long double C;
  long double E;
  long double F;
  long double B;
  free( D );
  long double DD;
  long double delta_x;
  long double delta_y;

  int count = 0;
  while ( count < 1000 && a > 0.001 && ( ( a - prev_max_delta ) / prev_max_delta >= tolerance || ( a - prev_max_delta ) / prev_max_delta <= -tolerance ) )
  {
    count++;
    prev_max_delta = a;
    for ( i = 0; i < N; i++ )
    {
      contributionx_m[i] = 0;
      contributiony_m[i] = 0;
    }
    for ( i = 1; i <= N; i++ )
    {
      if ( i == m )
      {
        continue;
      }
      contributionx_m[i - 1] = K[N * ( m - 1 ) + i - 1] * ( ( x[i - 1] - x[m - 1] ) - (long double)L[N * ( m - 1 ) + i - 1] * ( x[i - 1] - x[m - 1] ) / sqrtl( T[N * ( i - 1 ) + m - 1] ) );
      contributiony_m[i - 1] = K[N * ( m - 1 ) + i - 1] * ( ( y[i - 1] - y[m - 1] ) - (long double)L[N * ( m - 1 ) + i - 1] * ( y[i - 1] - y[m - 1] ) / sqrtl( T[N * ( i - 1 ) + m - 1] ) );
    }

    first_time = 1;
    long double prev_inner_max_delta = 0;
    while ( first_time || ( Delta[m - 1] > 0.000001 && ( prev_inner_max_delta - Delta[m - 1] ) / prev_inner_max_delta >= tolerance ) )
    {
      prev_inner_max_delta = Delta[m - 1];
      first_time = 0;
      C = -par_x[m - 1];
      E = -par_y[m - 1];

      F = 0;
      B = 0;
      DD = 0;
      for ( i = 1; i <= N; i++ )
      {
        if ( i == m )
        {
          continue;
        }
        F = F + K[N * ( i - 1 ) + m - 1] * ( 1 - (long double)L[N * ( i - 1 ) + m - 1] * powl( ( y[m - 1] - y[i - 1] ), 2 ) / powl( ( T[N * ( i - 1 ) + m - 1] ), 3.0 / 2 ) );
        B = B + K[N * ( i - 1 ) + m - 1] * ( (long double)L[N * ( i - 1 ) + m - 1] * ( x[m - 1] - x[i - 1] ) * ( y[m - 1] - y[i - 1] ) / powl( ( T[N * ( i - 1 ) + m - 1] ), 3.0 / 2 ) );
        DD = DD + K[N * ( i - 1 ) + m - 1] * ( 1 - (long double)L[N * ( i - 1 ) + m - 1] * powl( ( x[m - 1] - x[i - 1] ), 2 ) / powl( ( T[N * ( i - 1 ) + m - 1] ), 3.0 / 2 ) );
      }
      long double d1 = ( C * DD - E * B );
      long double d2 = ( F * DD - B * B );

      delta_x = (long double)d1 / d2;
      d1 = ( E * F - B * C );
      d2 = ( F * DD - B * B );
      delta_y = (long double)d1 / d2;

      x[m - 1] = x[m - 1] + delta_x;
      y[m - 1] = y[m - 1] + delta_y;
      for ( i = 1; i <= N; i++ )
      {
        if ( i == m )
        {
          continue;
        }
        T[N * ( m - 1 ) + i - 1] = powl( ( x[m - 1] - x[i - 1] ), 2 ) + powl( ( y[m - 1] - y[i - 1] ), 2 );
        T[N * ( i - 1 ) + m - 1] = T[N * ( m - 1 ) + i - 1];
      }

      par_x[m - 1] = 0;
      par_y[m - 1] = 0;
      for ( i = 1; i <= N; i++ )
      {
        if ( i == m )
        {
          continue;
        }
        par_x[m - 1] = par_x[m - 1] + K[N * ( i - 1 ) + m - 1] * ( ( x[m - 1] - x[i - 1] ) - (long double)L[N * ( i - 1 ) + m - 1] * ( x[m - 1] - x[i - 1] ) / sqrtl( T[N * ( i - 1 ) + m - 1] ) );
        par_y[m - 1] = par_y[m - 1] + K[N * ( i - 1 ) + m - 1] * ( ( y[m - 1] - y[i - 1] ) - (long double)L[N * ( i - 1 ) + m - 1] * ( y[m - 1] - y[i - 1] ) / sqrtl( T[N * ( i - 1 ) + m - 1] ) );
      }
      Delta[m - 1] = sqrtl( powl( ( par_x[m - 1] ), 2 ) + powl( ( par_y[m - 1] ), 2 ) );
    }
    for ( i = 1; i <= N; i++ )
    {
      if ( i == m )
      {
        continue;
      }
      par_x[i - 1] = par_x[i - 1] + K[N * ( m - 1 ) + i - 1] * ( ( x[i - 1] - x[m - 1] ) - (long double)L[N * ( m - 1 ) + i - 1] * ( x[i - 1] - x[m - 1] ) / sqrtl( T[N * ( i - 1 ) + m - 1] ) ) - contributionx_m[i - 1];
      par_y[i - 1] = par_y[i - 1] + K[N * ( m - 1 ) + i - 1] * ( ( y[i - 1] - y[m - 1] ) - (long double)L[N * ( m - 1 ) + i - 1] * ( y[i - 1] - y[m - 1] ) / sqrtl( T[N * ( i - 1 ) + m - 1] ) ) - contributiony_m[i - 1];
      Delta[i - 1] = sqrtl( powl( ( par_x[i - 1] ), 2 ) + powl( ( par_y[i - 1] ), 2 ) );
    }
    a = -1000000;
    m = -1;
    for ( i = 0; i < N; i++ )
    {
      if ( Delta[i] > a )
      {
        a = Delta[i];
        m = i + 1;
      }
    }
  }

  long double minx = 10000000;
  long double miny = 10000000;
  for ( i = 1; i <= N; i++ )
  {
    if ( x[i - 1] < minx )
    {
      minx = x[i - 1];
    }
    if ( y[i - 1] < miny )
    {
      miny = y[i - 1];
    }
  }

  for ( i = 1; i <= N; i++ )
  {
    x[i - 1] = x[i - 1] - minx;
    y[i - 1] = y[i - 1] - miny;
  }
  free( contributionx_m );
  free( contributiony_m );
  free( Delta );
  free( par_x );
  free( par_y );
  free( T );
  free( L );
  free( K );
  free( stack );
  free( checked );
  free( smap );
}

void Get_Spring_Data( int* Sixidx, int* deg, int* vadj, int* vidx, int N, int M, int** Clabel, int** Dlabel, int* Cidx, int Clabel_max, int maxdeg, int seed, int providecoords, long double* x, long double* y )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Get_Spring_Data                                                        */
  /*                                                                         */
  /*  If providecoords = 0, runs KKSL to get coordinates.                    */
  /*  If providecoords = 1, coordinates are provided.                        */
  /*                                                                         */
  /*  A straight-line embedding is created on the coordinates and a crossing */
  /*  list is worked out for the embedding.                                  */
  /*                                                                         */
  /***************************************************************************/

  double pi = 3.14159265358979323846;
  int i;
  int j;
  int k;
  double r;
  if ( providecoords == 0 )
  {
    KKSL( x, y, Sixidx, N, deg, vadj, vidx, 1, 10, seed );
    for ( i = 0; i < N; i++ )
    {
      r = (double)rand() / (double)RAND_MAX;
      x[i] = x[i] + 0.001 * r;
      r = (double)rand() / (double)RAND_MAX;
      y[i] = y[i] + 0.001 * r;
    }
  }

  int Clabel_count = 0;

  double s;
  double t;
  int v1;
  int v2;
  int v3;
  int v4;
  double a0;
  double a1;
  double a2;
  double a3;
  double b0;
  double b1;
  double b2;
  double b3;
  double xc;
  double yc;
  double ang_v1;
  double ang_v2;
  double ang_v3;
  int* Clabeltemp;
  int* Dlabeltemp;
  Clabeltemp = (int*)malloc( M * sizeof( int ) );
  Dlabeltemp = (int*)malloc( M * sizeof( int ) );

  int count;
  double* t_values;
  double* tsort;
  int* tcheck;
  t_values = (double*)malloc( M * sizeof( double ) );
  tsort = (double*)malloc( M * sizeof( double ) );
  tcheck = (int*)malloc( M * sizeof( int ) );

  for ( i = 1; i <= M; i++ )
  {
    count = 0;
    for ( j = 1; j <= M; j++ )
    {
      v1 = Sixidx[i - 1];
      v2 = Sixidx[M + i - 1];
      v3 = Sixidx[j - 1];
      v4 = Sixidx[M + j - 1];
      if ( v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4 )
      {
        continue;
      }
      a0 = x[v1 - 1];
      a1 = x[v2 - 1];
      a2 = x[v3 - 1];
      a3 = x[v4 - 1];
      b0 = y[v1 - 1];
      b1 = y[v2 - 1];
      b2 = y[v3 - 1];
      b3 = y[v4 - 1];
      if ( ( a1 - a0 < 0.0000001 && a0 - a1 < 0.0000001 ) && ( a2 - a3 < 0.0000001 && a3 - a2 < 0.0000001 ) )
      {
        // then parallel vertical.
        continue;
      }
      if ( ( ( b1 - b0 ) / ( a1 - a0 ) - ( b3 - b2 ) / ( a3 - a2 ) < 0.0000001 ) && ( ( b3 - b2 ) / ( a3 - a2 ) - ( b1 - b0 ) / ( a1 - a0 ) < 0.0000001 ) )
      {
        // then parallel.
        continue;
      }
      double t1 = ( ( a3 - a2 ) * ( b0 - b2 ) + ( b2 - b3 ) * ( a0 - a2 ) );
      double t2 = ( ( b3 - b2 ) * ( a1 - a0 ) + ( a2 - a3 ) * ( b1 - b0 ) );
      t = t1 / t2;
      if ( b3 - b2 >= 0.0000001 || b3 - b2 <= -0.0000001 )
      {
        s = ( (double)( b0 - b2 + ( b1 - b0 ) * t ) ) / ( (double)( b3 - b2 ) );
      }
      else
      {
        s = ( (double)( a0 - a2 + ( a1 - a0 ) * t ) ) / ( (double)( a3 - a2 ) );
      }
      if ( t > 0 && t < 1 && s > 0 && s < 1 )
      {
        count++;
        t_values[count - 1] = t;
        Clabel_count++;
        while ( Clabel_count >= Clabel_max )
        {
          *Clabel = (int*)realloc( *Clabel, ( Clabel_max + 100 ) * sizeof( int ) );
          *Dlabel = (int*)realloc( *Dlabel, ( Clabel_max + 100 ) * sizeof( int ) );
          Clabel_max = Clabel_max + 100;
        }
        ( *Clabel )[Clabel_count - 1] = j;

        xc = a0 + ( a1 - a0 ) * t;
        yc = b0 + ( b1 - b0 ) * t;
        ang_v1 = atan( ( b0 - yc ) / ( a0 - xc ) );
        if ( a0 < xc )
        {
          ang_v1 = ang_v1 + pi;
        }
        ang_v2 = atan( ( b1 - yc ) / ( a1 - xc ) );
        if ( a1 < xc )
        {
          ang_v2 = ang_v2 + pi;
        }
        ang_v3 = atan( ( b2 - yc ) / ( a2 - xc ) );
        if ( a2 < xc )
        {
          ang_v3 = ang_v3 + pi;
        }
        if ( ang_v2 < ang_v1 )
        {
          if ( ang_v3 > ang_v2 && ang_v3 < ang_v1 )
          {
            ( *Dlabel )[Clabel_count - 1] = -1;
          }
          else
          {
            ( *Dlabel )[Clabel_count - 1] = 1;
          }
        }
        else
        {
          if ( ang_v3 > ang_v2 || ang_v3 < ang_v1 )
          {
            ( *Dlabel )[Clabel_count - 1] = -1;
          }
          else
          {
            ( *Dlabel )[Clabel_count - 1] = 1;
          }
        }
      }
    }
    for ( j = 0; j < count; j++ )
    {
      tsort[j] = t_values[j];
      tcheck[j] = 0;
    }

    qsort( tsort, count, sizeof( double ), ascend_cmpfunc_double );
    int jc = 0;
    int tind = 0;
    for ( j = Clabel_count - count + 1; j <= Clabel_count; j++ )
    {
      jc++;
      tind = -1;
      for ( k = 0; k < count; k++ )
      {
        if ( t_values[k] == tsort[jc - 1] && tcheck[k] == 0 )
        {
          tcheck[k] = 1;
          tind = k + 1;
          break;
        }
      }
      if ( tind == -1 )
      {
        fprintf( stdout, "tind error \n" );
      }
      Clabeltemp[jc - 1] = ( *Clabel )[Clabel_count - count + tind - 1];
      Dlabeltemp[jc - 1] = ( *Dlabel )[Clabel_count - count + tind - 1];
    }
    jc = 0;
    for ( j = Clabel_count - count + 1; j <= Clabel_count; j++ )
    {
      jc++;
      ( *Clabel )[j - 1] = Clabeltemp[jc - 1];
      ( *Dlabel )[j - 1] = Dlabeltemp[jc - 1];
    }
    Cidx[i] = Clabel_count;
  }

  double* angles;
  int* vert_ind;
  double* vert_ang;
  angles = (double*)malloc( 2 * M * sizeof( double ) );
  vert_ind = (int*)malloc( N * N * sizeof( int ) );
  vert_ang = (double*)malloc( N * N * sizeof( double ) );

  double ang;
  double x1;
  double x2;
  double y1;
  double y2;
  int* next;
  int* prev;
  next = (int*)malloc( maxdeg * sizeof( int ) );
  prev = (int*)malloc( maxdeg * sizeof( int ) );
  int edge;

  for ( i = 1; i <= M; i++ )
  {
    v1 = Sixidx[i - 1];
    v2 = Sixidx[M + i - 1];
    x1 = x[v1 - 1];
    x2 = x[v2 - 1];
    y1 = y[v1 - 1];
    y2 = y[v2 - 1];
    if ( ( x2 - x1 ) < 0.0000001 && ( x1 - x2 ) < 0.0000001 )
    {
      if ( y1 > y2 )
      {
        angles[i - 1] = 3 * pi / 2;
        angles[M + i - 1] = pi / 2;
      }
      else
      {
        angles[i - 1] = pi / 2;
        angles[M + i - 1] = 3 * pi / 2;
      }
    }
    else
    {

      ang = atan( ( y2 - y1 ) / ( x2 - x1 ) );
      if ( x1 < x2 )
      {
        angles[i - 1] = ang;
        angles[M + i - 1] = ang + pi;
      }
      else
      {
        angles[i - 1] = ang + pi;
        angles[M + i - 1] = ang;
      }
    }
  }

  for ( i = 1; i <= N; i++ )
  {
    count = 0;
    for ( j = 1; j <= M; j++ )
    {
      if ( Sixidx[j - 1] == i )
      {
        count++;
        vert_ang[N * ( count - 1 ) + i - 1] = angles[j - 1];
        vert_ind[N * ( count - 1 ) + i - 1] = j;
      }
      if ( Sixidx[M + j - 1] == i )
      {
        count++;
        vert_ang[N * ( count - 1 ) + i - 1] = angles[M + j - 1];
        vert_ind[N * ( count - 1 ) + i - 1] = j;
      }
    }
    for ( j = 0; j < deg[i - 1]; j++ )
    {
      tsort[j] = vert_ang[N * ( j ) + i - 1];
      tcheck[j] = 0;
    }
    int jc = 0;
    int tind = 0;
    qsort( tsort, deg[i - 1], sizeof( double ), descend_cmpfunc_double );

    for ( j = 0; j < deg[i - 1]; j++ )
    {
      jc++;
      tind = -1;

      for ( k = 0; k < deg[i - 1]; k++ )
      {
        if ( vert_ang[N * ( k ) + i - 1] == tsort[jc - 1] && tcheck[k] == 0 )
        {
          tcheck[k] = 1;
          tind = k + 1;
          break;
        }
      }
      if ( tind == -1 )
      {
        fprintf( stdout, "tind error2 \n" );
      }
      Clabeltemp[jc - 1] = vert_ind[N * ( tind - 1 ) + i - 1];
    }
    for ( j = 0; j < deg[i - 1]; j++ )
    {
      vert_ind[N * ( j ) + i - 1] = Clabeltemp[j];
    }
  }
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= deg[i - 1]; j++ )
    {
      if ( j < deg[i - 1] )
      {
        next[j - 1] = j + 1;
      }
      else
      {
        next[j - 1] = 1;
      }
    }
    for ( j = 1; j <= deg[i - 1]; j++ )
    {
      if ( j == 1 )
      {
        prev[j - 1] = deg[i - 1];
      }
      else
      {
        prev[j - 1] = j - 1;
      }
    }
    for ( j = 1; j <= deg[i - 1]; j++ )
    {
      edge = vert_ind[N * ( j - 1 ) + i - 1];
      if ( Sixidx[edge - 1] == i )
      {
        Sixidx[2 * M + edge - 1] = vert_ind[N * ( next[j - 1] - 1 ) + i - 1];
        Sixidx[3 * M + edge - 1] = vert_ind[N * ( prev[j - 1] - 1 ) + i - 1];
      }
      else
      {
        Sixidx[4 * M + edge - 1] = vert_ind[N * ( next[j - 1] - 1 ) + i - 1];
        Sixidx[5 * M + edge - 1] = vert_ind[N * ( prev[j - 1] - 1 ) + i - 1];
      }
    }
  }

  free( Clabeltemp );
  free( Dlabeltemp );
  free( t_values );
  free( tsort );
  free( tcheck );
  free( angles );
  free( vert_ind );
  free( vert_ang );
  free( next );
  free( prev );
}

void UDFS( int v, int i, int* n, int* LP, int* Tcheck, int* T, int tcount, int* p, int* el, int* vidx, int* vadj, int* clst, int ccnt, int* eadj, int N, int* compidx, int compcnt, int* emap, int ecnt, int* elst, int* eidx, int* outputs, int M )
{

  /***************************************************************************/
  /*                                                                         */
  /*  UDFS                                                                   */
  /*                                                                         */
  /*  Upgraded DFS search for finding the maximal biconnected components of  */
  /*  a connected graph.                                                     */
  /*                                                                         */
  /***************************************************************************/

  n[v - 1] = i;
  LP[v - 1] = n[v - 1];
  i++;
  int w;
  int q;
  int e;
  int k;
  for ( w = vidx[v - 1] + 1; w <= vidx[v]; w++ )
  {
    q = vadj[w - 1];
    e = eadj[w - 1];
    if ( Tcheck[e - 1] == 0 )
    {
      Tcheck[e - 1] = 1;
      tcount++;
      T[tcount - 1] = e;
    }
    if ( n[q - 1] == 0 )
    {
      p[q - 1] = v;
      outputs[0] = tcount;
      outputs[1] = compcnt;
      outputs[2] = ccnt;
      outputs[3] = ecnt;

      UDFS( q, i, n, LP, Tcheck, T, tcount, p, el, vidx, vadj, clst, ccnt, eadj, N, compidx, compcnt, emap, ecnt, elst, eidx, outputs, M );
      tcount = outputs[0];
      compcnt = outputs[1];
      ccnt = outputs[2];
      ecnt = outputs[3];
      if ( LP[q - 1] >= n[v - 1] )
      { // then block completed
        if ( T[tcount - 1] != e )
        {
          int* vcheck;
          vcheck = (int*)calloc( N, sizeof( int ) );
          compcnt++;
          for ( k = tcount; k >= 1; k-- )
          {
            int v1 = el[emap[T[k - 1] - 1] - 1];
            int v2 = el[M + emap[T[k - 1] - 1] - 1];
            if ( vcheck[v1 - 1] == 0 )
            {
              vcheck[v1 - 1] = 1;
              ccnt++;
              clst[ccnt - 1] = v1;
            }
            if ( vcheck[v2 - 1] == 0 )
            {
              vcheck[v2 - 1] = 1;
              ccnt++;
              clst[ccnt - 1] = v2;
            }
            ecnt++;
            elst[ecnt - 1] = emap[T[k - 1] - 1];
            tcount = tcount - 1;
            if ( T[k - 1] == e )
            {
              break;
            }
          }
          free( vcheck );
          eidx[compcnt] = ecnt;
          compidx[compcnt] = ccnt;
        }
        else
        {
          tcount = tcount - 1;
        }
      }
      if ( LP[v - 1] > LP[q - 1] )
      {
        LP[v - 1] = LP[q - 1];
      }
    }
    else if ( q != p[v - 1] )
    {
      if ( LP[v - 1] > n[q - 1] )
      {
        LP[v - 1] = n[q - 1];
      }
    }
  }
  outputs[0] = tcount;
  outputs[1] = compcnt;
  outputs[2] = ccnt;
  outputs[3] = ecnt;
}

void UDFS5( int* n, int* LP, int i, int* Tcheck, int* T, int tcount, int* p, int* clst, int ccnt, int* compidx, int compcnt, int ecnt, int* elst, int* eidx, int* artverts, int M2, int* vadj2, int* vidx2, int* el2, int* eadj2, int v, int* el, int* vadj, int* vidx, int* eadj, int* outputs, int Mcomplimit, int artcnt, int N, int M )
{

  /***************************************************************************/
  /*                                                                         */
  /*  UDFS5                                                                  */
  /*                                                                         */
  /*  Upgraded DFS search for finding the maximal biconnected components of  */
  /*  a connected graph. Also creates a new biconnected graph from the       */
  /*  original graph. Designed for use with embed scheme 5.                  */
  /*                                                                         */
  /***************************************************************************/

  n[v - 1] = i;
  LP[v - 1] = n[v - 1];
  i++;
  int k;
  int u = 0;
  int w;
  int q;
  int e;
  int* vcheck;
  vcheck = (int*)calloc( N, sizeof( int ) );

  for ( w = vidx[v - 1] + 1; w <= vidx[v]; w++ )
  {
    q = vadj[w - 1];
    e = eadj[w - 1];
    if ( Tcheck[e - 1] == 0 )
    {
      Tcheck[e - 1] = 1;
      tcount++;
      T[tcount - 1] = e;
    }
    if ( n[q - 1] == 0 )
    {
      if ( u == 0 )
      {
        u = q;
      }
      p[q - 1] = v;
      int outputs2[8];
      UDFS5( n, LP, i, Tcheck, T, tcount, p, clst, ccnt, compidx, compcnt, ecnt, elst, eidx, artverts, M2, vadj2, vidx2, el2, eadj2, q, el, vadj, vidx, eadj, outputs2, Mcomplimit, artcnt, N, M );
      i = outputs2[0];
      tcount = outputs2[1];
      ccnt = outputs2[2];
      compcnt = outputs2[3];
      ecnt = outputs2[4];
      M2 = outputs2[5];
      Mcomplimit = outputs2[6];
      artcnt = outputs2[7];

      if ( LP[q - 1] >= n[v - 1] )
      {
        artcnt++;
        artverts[artcnt - 1] = v;
        if ( T[tcount - 1] != e )
        {
          for ( k = 1; k <= N; k++ )
          {
            vcheck[k - 1] == 0;
          }
          compcnt++;
          for ( k = tcount; k >= 1; k-- )
          {
            if ( vcheck[el[T[k - 1] - 1] - 1] == 0 )
            {

              vcheck[el[T[k - 1] - 1] - 1]++;
              ccnt++;
              clst[ccnt - 1] = el[T[k - 1] - 1];
            }
            if ( vcheck[el[M + T[k - 1] - 1] - 1] == 0 )
            {
              vcheck[el[M + T[k - 1] - 1] - 1]++;
              ccnt++;
              clst[ccnt - 1] = el[M + T[k - 1] - 1];
            }
            ecnt++;
            elst[ecnt - 1] = T[k - 1];
            tcount--;
            if ( T[k - 1] == e )
            {
              break;
            }
          }
          eidx[compcnt] = ecnt;
          compidx[compcnt] = ccnt;
        }
        else
        {
          tcount--;
        }
        if ( q == u && p[v - 1] > 0 )
        {
          M2++;
          if ( M2 > Mcomplimit )
          {
            // realloc the new lists.
            el2 = (int*)realloc( el2, ( M2 + 20 ) * sizeof( int ) );
            vadj2 = (int*)realloc( vadj2, 2 * ( M2 + 20 ) * sizeof( int ) );
            eadj2 = (int*)realloc( eadj2, 2 * ( M2 + 20 ) * sizeof( int ) );
            int rr;
            int s;
            for ( rr = 1; rr >= 0; rr-- )
            {
              for ( s = Mcomplimit; s >= 1; s-- )
              {
                el2[rr * ( M2 + 20 ) + s - 1] = el2[rr * Mcomplimit + s - 1];
              }
            }
            Mcomplimit = M2 + 20;
          }
          if ( q < p[v - 1] )
          {
            el2[M2 - 1] = q;
            el2[Mcomplimit + M2 - 1] = p[v - 1];
          }
          else
          {
            el2[M2 - 1] = p[v - 1];
            el2[Mcomplimit + M2 - 1] = q;
          }
          for ( k = vidx2[N]; k >= vidx2[q] + 1; k-- )
          {
            vadj2[k] = vadj2[k - 1];
            eadj2[k] = eadj2[k - 1];
          }
          vadj2[vidx2[q]] = p[v - 1];
          eadj2[vidx2[q]] = M2;
          for ( k = q + 1; k <= N + 1; k++ )
          {
            vidx2[k - 1]++;
          }
          for ( k = vidx2[N]; k >= vidx2[p[v - 1]] + 1; k-- )
          {
            vadj2[k] = vadj2[k - 1];
            eadj2[k] = eadj2[k - 1];
          }
          vadj2[vidx2[p[v - 1]]] = q;
          eadj2[vidx2[p[v - 1]]] = M2;
          for ( k = p[v - 1] + 1; k <= N + 1; k++ )
          {
            vidx2[k - 1]++;
          }
        }
        if ( q != u )
        {
          M2++;
          if ( M2 > Mcomplimit )
          {
            // realloc the new lists.
            el2 = (int*)realloc( el2, ( M2 + 20 ) * sizeof( int ) );
            vadj2 = (int*)realloc( vadj2, 2 * ( M2 + 20 ) * sizeof( int ) );
            eadj2 = (int*)realloc( eadj2, 2 * ( M2 + 20 ) * sizeof( int ) );
            int rr;
            int s;
            for ( rr = 1; rr >= 0; rr-- )
            {
              for ( s = Mcomplimit; s >= 1; s-- )
              {
                el2[rr * ( M2 + 20 ) + s - 1] = el2[rr * Mcomplimit + s - 1];
              }
            }
            Mcomplimit = M2 + 20;
          }
          if ( q < u )
          {
            el2[M2 - 1] = q;
            el2[Mcomplimit + M2 - 1] = u;
          }
          else
          {
            el2[M2 - 1] = u;
            el2[Mcomplimit + M2 - 1] = q;
          }
          for ( k = vidx2[N]; k >= vidx2[q] + 1; k-- )
          {
            vadj2[k] = vadj2[k - 1];
            eadj2[k] = eadj2[k - 1];
          }
          vadj2[vidx2[q]] = u;
          eadj2[vidx2[q]] = M2;
          for ( k = q + 1; k <= N + 1; k++ )
          {
            vidx2[k - 1]++;
          }
          for ( k = vidx2[N]; k >= vidx2[u] + 1; k-- )
          {
            vadj2[k] = vadj2[k - 1];
            eadj2[k] = eadj2[k - 1];
          }
          vadj2[vidx2[u]] = q;
          eadj2[vidx2[u]] = M2;
          for ( k = u + 1; k <= N + 1; k++ )
          {
            vidx2[k - 1]++;
          }
        }
      }
      if ( LP[v - 1] > LP[q - 1] )
      {
        LP[v - 1] = LP[q - 1];
      }
    }
    else if ( q != p[v - 1] )
    {
      if ( LP[v - 1] > n[q - 1] )
      {
        LP[v - 1] = n[q - 1];
      }
    }
  }
  free( vcheck );
  outputs[0] = i;
  outputs[1] = tcount;
  outputs[2] = ccnt;
  outputs[3] = compcnt;
  outputs[4] = ecnt;
  outputs[5] = M2;
  outputs[6] = Mcomplimit;
  outputs[7] = artcnt;
}

void UDFS4( int v, int i, int* n, int* LP, int* Tcheck, int* T, int tcount, int* p, int* el, int* vidx, int* vadj, int* clst, int ccnt, int* eadj, int N, int* compidx, int compcnt, int ecnt, int* elst, int* eidx, int* artvert, int artcnt, int* outputs, int M )
{

  /***************************************************************************/
  /*                                                                         */
  /*  UDFS4                                                                  */
  /*                                                                         */
  /*  Upgraded DFS search for finding the maximal biconnected components of  */
  /*  a connected graph. Does not create a modified biconnected graph at the */
  /*  same time. Designed for use with embed scheme 4.                       */
  /*                                                                         */
  /***************************************************************************/

  n[v - 1] = i;
  LP[v - 1] = n[v - 1];
  i++;
  int u = 0;
  int w;
  int q;
  int e;
  int k;

  for ( w = vidx[v - 1] + 1; w <= vidx[v]; w++ )
  {
    q = vadj[w - 1];
    e = eadj[w - 1];
    if ( Tcheck[e - 1] == 0 )
    {
      Tcheck[e - 1] = 1;
      tcount++;
      T[tcount - 1] = e;
    }
    if ( n[q - 1] == 0 )
    {
      if ( u == 0 )
      {
        u = q;
      }
      p[q - 1] = v;
      int outputs2[6];
      UDFS4( q, i, n, LP, Tcheck, T, tcount, p, el, vidx, vadj, clst, ccnt, eadj, N, compidx, compcnt, ecnt, elst, eidx, artvert, artcnt, outputs2, M );
      i = outputs2[0];
      tcount = outputs2[1];
      ccnt = outputs2[2];
      compcnt = outputs2[3];
      ecnt = outputs2[4];
      artcnt = outputs2[5];

      if ( LP[q - 1] >= n[v - 1] )
      {
        artcnt++;
        artvert[artcnt - 1] = v;
        if ( T[tcount - 1] != e )
        {
          int* vcheck;
          vcheck = (int*)calloc( N, sizeof( int ) );
          compcnt++;
          for ( k = tcount; k >= 1; k-- )
          {
            if ( vcheck[el[T[k - 1] - 1] - 1] == 0 )
            {
              vcheck[el[T[k - 1] - 1] - 1]++;
              ccnt++;
              clst[ccnt - 1] = el[T[k - 1] - 1];
            }
            if ( vcheck[el[M + T[k - 1] - 1] - 1] == 0 )
            {
              vcheck[el[M + T[k - 1] - 1] - 1]++;
              ccnt++;
              clst[ccnt - 1] = el[M + T[k - 1] - 1];
            }
            ecnt++;
            elst[ecnt - 1] = T[k - 1];
            tcount--;
            if ( T[k - 1] == e )
            {
              break;
            }
          }
          eidx[compcnt] = ecnt;
          compidx[compcnt] = ccnt;
          free( vcheck );
        }
        else
        {
          tcount--;
        }
      }
      if ( LP[v - 1] > LP[q - 1] )
      {
        LP[v - 1] = LP[q - 1];
      }
    }
    else if ( q != p[v - 1] )
    {
      if ( LP[v - 1] > n[q - 1] )
      {
        LP[v - 1] = n[q - 1];
      }
    }
  }
  outputs[0] = i;
  outputs[1] = tcount;
  outputs[2] = ccnt;
  outputs[3] = compcnt;
  outputs[4] = ecnt;
  outputs[5] = artcnt;
}

void DFSrelabel4( int* n, int* LP, int i, int* p, int vcnt, int* LP2, int* backedge, int* eadj, int v, int* vidx, int* vadj, int N, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  DFSrelabel4                                                            */
  /*                                                                         */
  /*  Used to relabel vertices for the DFS. Used for embed scheme 4.         */
  /*                                                                         */
  /***************************************************************************/

  int q;
  int w;
  n[v - 1] = i;
  LP[v - 1] = n[v - 1];
  LP2[v - 1] = n[v - 1];
  i++;
  int outputs2[2];
  for ( w = vidx[v - 1] + 1; w <= vidx[v]; w++ )
  {
    q = vadj[w - 1];
    if ( n[q - 1] == 0 )
    {
      p[q - 1] = v;
      vcnt++;
      int outputs2[2];
      DFSrelabel4( n, LP, i, p, vcnt, LP2, backedge, eadj, q, vidx, vadj, N, outputs2 );
      i = outputs2[0];
      vcnt = outputs2[1];
      if ( LP[v - 1] > LP[q - 1] )
      {
        LP[v - 1] = LP[q - 1];
      }
    }
    else if ( q != p[v - 1] )
    {
      if ( LP[v - 1] > n[q - 1] )
      {
        LP[v - 1] = n[q - 1];
      }
      backedge[eadj[w - 1] - 1] = 1;
    }
  }
  for ( w = vidx[v - 1] + 1; w <= vidx[v]; w++ )
  {
    q = vadj[w - 1];
    if ( p[q - 1] == v )
    {
      if ( LP[q - 1] != LP[v - 1] && LP[q - 1] < LP2[v - 1] )
      {
        LP2[v - 1] = LP[q - 1];
      }
      if ( LP2[q - 1] < LP2[v - 1] )
      {
        LP2[v - 1] = LP2[q - 1];
      }
    }
    else
    {
      if ( LP[v - 1] != n[q - 1] && n[q - 1] < LP2[v - 1] )
      {
        LP2[v - 1] = n[q - 1];
      }
    }
  }
  outputs[0] = i;
  outputs[1] = vcnt;
}

void Segment_Bipartite( int h, int wj, int M, int loz, int* Alefts, int* alfirst, int* alnext, int* Arights, int* arfirst, int* arnext, int* Elefts, int* elfirst, int* elnext, int* Erights, int* erfirst, int* ernext, int* Ae, int Aecnt, int* Ee, int N, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Segment_Bipartite                                                      */
  /*                                                                         */
  /*  Confirms a segment is bipartite. Used in the planarity algorithm.      */
  /*                                                                         */
  /***************************************************************************/

  // the outputs vector contains d, inside_check, planar.
  int r;
  int s;
  int planar = 1;
  int ell = h + 1;
  int mcnt = 0;
  int mcnt2 = 0;
  int maxarb = 0;
  int maxalb = 0;
  int qv;
  r = 0;

  while ( r < Aecnt )
  {
    r++;
    if ( Ae[r - 1] == wj )
    {
      for ( s = r; s <= Aecnt - 1; s++ )
      {
        Ae[s - 1] = Ae[s];
        Ee[s - 1] = Ee[s];
      }
      Aecnt--;
      r--;
    }
  }
  int i = h;
  planar = 1;
  if ( i == 0 )
  {
    alfirst[0] = 1;
    elfirst[0] = 1;
    for ( r = 1; r <= Aecnt; r++ )
    {
      Alefts[r - 1] = Ae[r - 1];
      Elefts[r - 1] = Ee[r - 1];
      if ( r < Aecnt )
      {
        alnext[r - 1] = r + 1;
        elnext[r - 1] = r + 1;
      }
    }
    i = 1;
    outputs[0] = 1;
    outputs[1] = 1;
    outputs[5] = Aecnt;
    outputs[6] = 0;
    return;
  }

  int* aL = (int*)calloc( M, sizeof( int ) );
  int* eL = (int*)calloc( M, sizeof( int ) );
  int* aR = (int*)calloc( M, sizeof( int ) );
  int* eR = (int*)calloc( M, sizeof( int ) );
  int alcnt = Aecnt;
  int arcnt = 0;
  int elcnt = Aecnt;
  int ercnt = 0;
  for ( r = 1; r <= Aecnt; r++ )
  {
    aL[r - 1] = Ae[r - 1];
    eL[r - 1] = Ee[r - 1];
  }

  qv = alfirst[i - 1];
  for ( r = 1; r <= N; r++ )
  {
    if ( qv == 0 )
    {
      break;
    }
    if ( Alefts[qv - 1] > maxalb )
    {
      maxalb = Alefts[qv - 1];
    }
    qv = alnext[qv - 1];
  }
  qv = arfirst[i - 1];
  for ( r = 1; r <= N; r++ )
  {
    if ( qv == 0 )
    {
      break;
    }
    if ( Arights[qv - 1] > maxarb )
    {
      maxarb = Arights[qv - 1];
    }
    qv = arnext[qv - 1];
  }

  while ( i > 0 && ( maxalb > loz || maxarb > loz ) )
  {
    if ( alfirst[i - 1] > 0 && maxalb > loz )
    {
      if ( arfirst[i - 1] > 0 && maxarb > loz )
      {
        planar = 0;
        outputs[1] = planar;
        free( aL );
        free( aR );
        free( eL );
        free( eR );
        return;
      }
      qv = alfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        arcnt++;
        aR[arcnt - 1] = Alefts[qv - 1];
        qv = alnext[qv - 1];
      }
      qv = arfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        alcnt++;
        aL[alcnt - 1] = Arights[qv - 1];
        qv = arnext[qv - 1];
      }
      qv = elfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        ercnt++;
        eR[ercnt - 1] = Elefts[qv - 1];
        qv = elnext[qv - 1];
      }
      qv = erfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        elcnt++;
        eL[elcnt - 1] = Erights[qv - 1];
        qv = ernext[qv - 1];
      }
    }
    else
    {
      qv = alfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        alcnt++;
        aL[alcnt - 1] = Alefts[qv - 1];
        qv = alnext[qv - 1];
      }
      qv = arfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        arcnt++;
        aR[arcnt - 1] = Arights[qv - 1];
        qv = arnext[qv - 1];
      }
      qv = elfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        elcnt++;
        eL[elcnt - 1] = Elefts[qv - 1];
        qv = elnext[qv - 1];
      }
      qv = erfirst[i - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        ercnt++;
        eR[ercnt - 1] = Erights[qv - 1];
        qv = ernext[qv - 1];
      }
    }

    i--;
    if ( i == 0 )
    {
      break;
    }
    maxarb = 0;
    maxalb = 0;

    qv = alfirst[i - 1];
    for ( r = 1; r <= N; r++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( Alefts[qv - 1] > maxalb )
      {
        maxalb = Alefts[qv - 1];
      }
      qv = alnext[qv - 1];
    }
    qv = arfirst[i - 1];
    for ( r = 1; r <= N; r++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( Arights[qv - 1] > maxarb )
      {
        maxarb = Arights[qv - 1];
      }
      qv = arnext[qv - 1];
    }
  }

  int* Alefts2;
  int* Arights2;
  int* Elefts2;
  int* Erights2;
  int* alfirst2;
  int* arfirst2;
  int* erfirst2;
  int* elfirst2;
  int* alnext2;
  int* arnext2;
  int* ernext2;
  int* elnext2;
  int free2arb = 1;
  int free2alb = 1;
  int free2elb = 1;
  int free2erb = 1;
  Alefts2 = (int*)calloc( M, sizeof( int ) );
  Arights2 = (int*)calloc( M, sizeof( int ) );
  Elefts2 = (int*)calloc( M, sizeof( int ) );
  Erights2 = (int*)calloc( M, sizeof( int ) );
  alfirst2 = (int*)calloc( M, sizeof( int ) );
  arfirst2 = (int*)calloc( M, sizeof( int ) );
  erfirst2 = (int*)calloc( M, sizeof( int ) );
  elfirst2 = (int*)calloc( M, sizeof( int ) );
  alnext2 = (int*)calloc( M, sizeof( int ) );
  arnext2 = (int*)calloc( M, sizeof( int ) );
  ernext2 = (int*)calloc( M, sizeof( int ) );
  elnext2 = (int*)calloc( M, sizeof( int ) );
  int prevfree2 = 0;

  for ( r = 1; r <= i; r++ )
  {
    qv = alfirst[r - 1]; // lefts..
    prevfree2 = -1;
    for ( s = 1; s <= N; s++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( alfirst2[r - 1] == 0 )
      {
        alfirst2[r - 1] = free2alb;
      }
      else
      {
        alnext2[prevfree2 - 1] = free2alb;
      }
      Alefts2[free2alb - 1] = Alefts[qv - 1];
      prevfree2 = free2alb;
      free2alb++;
      qv = alnext[qv - 1];
    }
    qv = elfirst[r - 1]; // left edges.
    prevfree2 = -1;

    for ( s = 1; s <= N; s++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( elfirst2[r - 1] == 0 )
      {
        elfirst2[r - 1] = free2elb;
      }
      else
      {
        elnext2[prevfree2 - 1] = free2elb;
      }
      Elefts2[free2elb - 1] = Elefts[qv - 1];
      prevfree2 = free2elb;
      free2elb++;
      qv = elnext[qv - 1];
    }
    qv = arfirst[r - 1]; // right..
    prevfree2 = -1;
    for ( s = 1; s <= N; s++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( arfirst2[r - 1] == 0 )
      {
        arfirst2[r - 1] = free2arb;
      }
      else
      {
        arnext2[prevfree2 - 1] = free2arb;
      }
      Arights2[free2arb - 1] = Arights[qv - 1];
      prevfree2 = free2arb;
      free2arb++;
      qv = arnext[qv - 1];
    }
    qv = erfirst[r - 1]; // right edges.
    prevfree2 = -1;
    for ( s = 1; s <= N; s++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( erfirst2[r - 1] == 0 )
      {
        erfirst2[r - 1] = free2erb;
      }
      else
      {
        ernext2[prevfree2 - 1] = free2erb;
      }
      Erights2[free2erb - 1] = Erights[qv - 1];
      prevfree2 = free2erb;
      free2erb++;
      qv = ernext[qv - 1];
    }
  }
  i++;

  if ( alcnt > 0 )
  {
    alfirst2[i - 1] = free2alb; // lefts..
  }
  for ( s = 1; s <= alcnt; s++ )
  {
    Alefts2[free2alb - 1] = aL[s - 1];
    if ( s < alcnt )
    {
      alnext2[free2alb - 1] = free2alb + 1;
    }
    free2alb++;
  }
  if ( elcnt > 0 )
  {
    elfirst2[i - 1] = free2elb; // lefts..
  }
  for ( s = 1; s <= elcnt; s++ )
  {
    Elefts2[free2elb - 1] = eL[s - 1];
    if ( s < elcnt )
    {
      elnext2[free2elb - 1] = free2elb + 1;
    }
    free2elb++;
  }
  if ( arcnt > 0 )
  {
    arfirst2[i - 1] = free2arb; // rights..
  }
  for ( s = 1; s <= arcnt; s++ )
  {
    Arights2[free2arb - 1] = aR[s - 1];
    if ( s < arcnt )
    {
      arnext2[free2arb - 1] = free2arb + 1;
    }
    free2arb++;
  }
  if ( ercnt > 0 )
  {
    erfirst2[i - 1] = free2erb; // rights..
  }
  // fprintf(stdout,"eR %d i %d \n",eR[0],i);
  for ( s = 1; s <= ercnt; s++ )
  {
    Erights2[free2erb - 1] = eR[s - 1];
    if ( s < ercnt )
    {
      ernext2[free2erb - 1] = free2erb + 1;
    }
    free2erb++;
  }

  // copy back to orig albs and arbs.
  for ( r = 1; r <= M; r++ )
  {
    Alefts[r - 1] = Alefts2[r - 1];
    Arights[r - 1] = Arights2[r - 1];
    Elefts[r - 1] = Elefts2[r - 1];
    Erights[r - 1] = Erights2[r - 1];
    arfirst[r - 1] = arfirst2[r - 1];
    alfirst[r - 1] = alfirst2[r - 1];
    erfirst[r - 1] = erfirst2[r - 1];
    elfirst[r - 1] = elfirst2[r - 1];
    arnext[r - 1] = arnext2[r - 1];
    alnext[r - 1] = alnext2[r - 1];
    ernext[r - 1] = ernext2[r - 1];
    elnext[r - 1] = elnext2[r - 1];
    // nextfreealb[r-1]=r;
    // nextfreearb[r-1]=r;
  }

  outputs[0] = i;
  outputs[1] = planar;
  outputs[3] = free2elb;
  outputs[4] = free2erb;
  outputs[5] = free2arb;
  outputs[6] = free2alb;

  free( aL );
  free( aR );
  free( eL );
  free( eR );
  free( Alefts2 );
  free( Arights2 );
  free( Elefts2 );
  free( Erights2 );
  free( alfirst2 );
  free( arfirst2 );
  free( erfirst2 );
  free( elfirst2 );
  free( alnext2 );
  free( arnext2 );
  free( ernext2 );
  free( elnext2 );
}

void Strongly_Planar( int* L, int* EL, int* backedge, int e0, int* p2, int* neweadj, int* newvadj, int* newvidx, int* newel, int M, int N, int* newdeg, int* LP2, int* alpha, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Strongly_Planar                                                        */
  /*                                                                         */
  /*  Checks if a subgraph is strongly planar. Used in planarity checker.    */
  /*                                                                         */
  /***************************************************************************/

  // outputs are Lstartindex,planar
  int r;
  int i;
  int j;
  int wj;
  int wjminus;
  int ej;
  int lcnt = 0;
  int elcnt = 0;
  int planar = 1;
  int x;
  int y;
  if ( backedge[e0 - 1] == 0 )
  {
    x = newel[e0 - 1];
    y = newel[M + e0 - 1];
  }
  else
  {
    x = newel[M + e0 - 1];
    y = newel[e0 - 1];
  }
  int* Spinepath;
  Spinepath = (int*)calloc( N, sizeof( int ) );
  int currv = y;
  int curre = e0;
  int scnt = 0;
  int w0;
  if ( backedge[e0 - 1] == 0 )
  {
    for ( i = 1; i <= M; i++ )
    {
      if ( backedge[curre - 1] == 1 )
      {
        w0 = currv;
        break;
      }
      scnt++;

      Spinepath[scnt - 1] = currv;
      curre = neweadj[newvidx[currv - 1]];
      currv = newvadj[newvidx[currv - 1]];
    }
  }
  else
  {
    curre = e0;
    w0 = y;
  }

  currv = x;
  int w1 = y;
  for ( i = 1; i <= M; i++ )
  {
    if ( currv == w0 )
    {
      break;
    }
    w1 = currv;
    currv = p2[currv - 1];
  }

  int* B;
  int h = 0;
  // all linked lists now
  int* Alefts;
  int* Arights;
  int* Elefts;
  int* Erights;
  int* alfirst;
  int* alnext;
  int* arfirst;
  int* arnext;
  int* elfirst;
  int* elnext;
  int* erfirst;
  int* ernext;
  int freear = 0;
  int freeal = 0;
  int freeer = 0;
  int freeel = 0;

  B = (int*)calloc( M, sizeof( int ) );
  Alefts = (int*)calloc( M, sizeof( int ) );
  Arights = (int*)calloc( M, sizeof( int ) );
  Elefts = (int*)calloc( M, sizeof( int ) );
  Erights = (int*)calloc( M, sizeof( int ) );
  alfirst = (int*)calloc( M, sizeof( int ) );
  alnext = (int*)calloc( M, sizeof( int ) );
  arfirst = (int*)calloc( M, sizeof( int ) );
  arnext = (int*)calloc( M, sizeof( int ) );
  elfirst = (int*)calloc( M, sizeof( int ) );
  elnext = (int*)calloc( M, sizeof( int ) );
  erfirst = (int*)calloc( M, sizeof( int ) );
  ernext = (int*)calloc( M, sizeof( int ) );
  // int *nextfreear;
  // int *nextfreeal;
  // int nextfreecntarb=0;
  // int nextfreecntalb=0;
  // nextfreear = (int*)calloc(M,sizeof(int));
  // nextfreeal = (int*)calloc(M,sizeof(int));
  // for(r=1; r<=M; r++){
  // nextfreear[r-1]=r;
  // nextfreeal[r-1]=r;
  //}
  int* Ae;
  Ae = (int*)calloc( N, sizeof( int ) );
  int Aecnt = 0;
  int* Ee;
  Ee = (int*)calloc( N, sizeof( int ) );
  int Eecnt = 0;
  int loz;
  int maxmalbarb;
  int outputs2[7];
  int d;
  int inside_check;
  int qv;
  int prevqv;

  L[0] = w0;
  lcnt = 1;

  for ( j = scnt; j >= 1; j-- )
  {
    wj = Spinepath[j - 1];
    if ( j > 1 )
    {
      wjminus = Spinepath[j - 2];
    }
    else
    {
      wjminus = p2[wj - 1];
    }
    for ( r = 2; r <= newdeg[wj - 1]; r++ )
    {
      ej = neweadj[newvidx[wj - 1] + r - 1];
      // fprintf(stdout,"deg: %d wj: %d\n",newdeg[wj-1],wj);
      Strongly_Planar( Ae, Ee, backedge, ej, p2, neweadj, newvadj, newvidx, newel, M, N, newdeg, LP2, alpha, outputs2 );
      // fprintf(stdout,"r1 \n");
      Aecnt = outputs2[0];
      Eecnt = outputs2[1];
      planar = outputs2[2];
      if ( planar == 0 )
      {
        // fprintf(stdout,"exiting non-planar1 \n");
        outputs[2] = 0;
        free( Spinepath );
        free( Alefts );
        free( Arights );
        free( Elefts );
        free( Erights );
        free( alfirst );
        free( arfirst );
        free( erfirst );
        free( elfirst );
        free( arnext );
        free( alnext );
        free( ernext );
        free( elnext );
        // free(nextfreear);
        // free(nextfreeal);
        free( Ae );
        free( Ee );
        return;
      }
      if ( backedge[ej - 1] == 1 )
      {
        loz = newvadj[newvidx[wj - 1] + r - 1];
      }
      else
      {
        loz = LP2[newvadj[newvidx[wj - 1] + r - 1] - 1];
      }

      Segment_Bipartite( h, wj, M, loz, Alefts, alfirst, alnext, Arights, arfirst, arnext, Elefts, elfirst, elnext, Erights, erfirst, ernext, Ae, Aecnt, Ee, N, outputs2 );
      h = outputs2[0];
      planar = outputs2[1];

      if ( planar == 0 )
      {
        outputs[2] = 0;
        // fprintf(stdout,"exiting non-planar2 \n");
        return;
      }
    }

    if ( h > 0 )
    {
      int maxalbarb = 0;
      qv = alfirst[h - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        if ( Alefts[qv - 1] > maxalbarb )
        {
          maxalbarb = Alefts[qv - 1];
        }
        qv = alnext[qv - 1];
      }
      qv = arfirst[h - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        if ( Arights[qv - 1] > maxalbarb )
        {
          maxalbarb = Arights[qv - 1];
        }
        qv = arnext[qv - 1];
      }
      // fprintf(stdout,"check h: %d maxalbarb: %d wjminus: %d \n",h,maxalbarb,wjminus);
      while ( ( wjminus > 0 ) && h > 0 && ( alfirst[h - 1] > 0 || arfirst[h - 1] > 0 ) && ( maxalbarb == wjminus ) )
      {
        qv = alfirst[h - 1];
        for ( r = 1; r <= N; r++ )
        {
          if ( qv == 0 )
          {
            break;
          }
          if ( Alefts[qv - 1] == wjminus )
          {
            // nextfreecntalb--;
            // fprintf(stdout,"nxfc %d \n",nextfreecntalb);
            // nextfreeal[nextfreecntalb]=qv;
            if ( alfirst[h - 1] != qv )
            {
              alnext[prevqv - 1] = alnext[qv - 1];
            }
            else
            {
              alfirst[h - 1] = alnext[qv - 1];
            }
            alnext[qv - 1] = 0;
            Alefts[qv - 1] = 0;
          }
          prevqv = qv;
          qv = alnext[qv - 1];
        }
        qv = arfirst[h - 1];
        for ( r = 1; r <= N; r++ )
        {
          if ( qv == 0 )
          {
            break;
          }
          if ( Arights[qv - 1] == wjminus )
          {
            // nextfreecntarb--;
            // nextfreear[nextfreecntarb]=qv;
            if ( arfirst[h - 1] != qv )
            {
              arnext[prevqv - 1] = arnext[qv - 1];
            }
            else
            {
              arfirst[h - 1] = arnext[qv - 1];
            }
            arnext[qv - 1] = 0;
            Arights[qv - 1] = 0;
          }
          prevqv = qv;
          qv = arnext[qv - 1];
        }
        if ( alfirst[h - 1] == 0 && arfirst[h - 1] == 0 )
        {
          int qvold;
          qv = elfirst[h - 1];
          elfirst[h - 1] = 0;
          for ( r = 1; r <= N; r++ )
          {
            if ( qv == 0 )
            {
              break;
            }
            // fprintf(stdout,"k1: %d \n",Elefts[qv-1]);
            alpha[Elefts[qv - 1] - 1] = 1;
            qvold = qv;
            qv = elnext[qv - 1];
            elnext[qvold - 1] = 0;
          }
          qv = erfirst[h - 1];
          erfirst[h - 1] = 0;
          for ( r = 1; r <= N; r++ )
          {
            if ( qv == 0 )
            {
              break;
            }
            // fprintf(stdout,"-k2: %d \n",Erights[qv-1]);
            alpha[Erights[qv - 1] - 1] = -1;
            qvold = qv;
            qv = ernext[qv - 1];
            ernext[qvold - 1] = 0;
          }
          // remove the block of elefts and erights

          h--;
          if ( h == 0 )
          {
            break;
          }
        }
        // find new maxalarb
        maxalbarb = 0;
        qv = alfirst[h - 1];
        for ( r = 1; r <= N; r++ )
        {
          if ( qv == 0 )
          {
            break;
          }
          if ( Alefts[qv - 1] > maxalbarb )
          {
            maxalbarb = Alefts[qv - 1];
          }
          qv = alnext[qv - 1];
        }
        qv = arfirst[h - 1];
        for ( r = 1; r <= N; r++ )
        {
          if ( qv == 0 )
          {
            break;
          }
          if ( Arights[qv - 1] > maxalbarb )
          {
            maxalbarb = Arights[qv - 1];
          }
          qv = arnext[qv - 1];
        }
      }
    }
  }

  int maxarb;
  int maxalb;
  int ell;

  for ( ell = 1; ell <= h; ell++ )
  {
    maxalb = 0;
    maxarb = 0;
    qv = alfirst[ell - 1];
    for ( r = 1; r <= N; r++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( Alefts[qv - 1] > maxalb )
      {
        maxalb = Alefts[qv - 1];
      }
      qv = alnext[qv - 1];
    }
    qv = arfirst[ell - 1];
    for ( r = 1; r <= N; r++ )
    {
      if ( qv == 0 )
      {
        break;
      }
      if ( Arights[qv - 1] > maxarb )
      {
        maxarb = Arights[qv - 1];
      }
      qv = arnext[qv - 1];
    }
    if ( ( alfirst[ell - 1] > 0 && maxalb >= w1 ) && ( arfirst[ell - 1] > 0 && maxarb >= w1 ) )
    {
      free( Spinepath );
      free( Alefts );
      free( Arights );
      free( Elefts );
      free( Erights );
      free( alfirst );
      free( arfirst );
      free( erfirst );
      free( elfirst );
      free( arnext );
      free( alnext );
      free( ernext );
      free( elnext );
      // free(nextfreear);
      // free(nextfreeal);
      free( Ae );
      free( Ee );
      outputs[2] = 0;
      return;
    }
    if ( maxarb >= w1 )
    {
      qv = elfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        // fprintf(stdout,"-k3: %d \n",Elefts[qv-1]);
        alpha[Elefts[qv - 1] - 1] = -1;
        qv = elnext[qv - 1];
      }
      qv = erfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        // fprintf(stdout,"k4: %d \n",Erights[qv-1]);
        alpha[Erights[qv - 1] - 1] = 1;
        qv = ernext[qv - 1];
      }
    }
    else
    {
      qv = elfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        // fprintf(stdout,"k5: %d \n",Elefts[qv-1]);
        alpha[Elefts[qv - 1] - 1] = 1;
        qv = elnext[qv - 1];
      }
      qv = erfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        // fprintf(stdout,"-k6: %d \n",Erights[qv-1]);
        alpha[Erights[qv - 1] - 1] = -1;
        qv = ernext[qv - 1];
      }
    }
    if ( alfirst[ell - 1] > 0 && maxalb >= w1 )
    {
      // fprintf(stdout,"loop1 \n");
      qv = arfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        lcnt++;
        L[lcnt - 1] = Arights[qv - 1];
        qv = arnext[qv - 1];
      }
      qv = alfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        lcnt++;
        L[lcnt - 1] = Alefts[qv - 1];
        qv = alnext[qv - 1];
      }
    }
    else
    {
      // fprintf(stdout,"loop2 \n");
      qv = alfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        lcnt++;
        L[lcnt - 1] = Alefts[qv - 1];
        qv = alnext[qv - 1];
      }
      qv = arfirst[ell - 1];
      for ( r = 1; r <= N; r++ )
      {
        if ( qv == 0 )
        {
          break;
        }
        lcnt++;
        L[lcnt - 1] = Arights[qv - 1];
        qv = arnext[qv - 1];
      }
    }
  }

  for ( i = 1; i <= lcnt; i++ )
  {
    EL[i - 1] = e0;
  }

  // maybe L need sorted here.
  outputs[0] = lcnt;
  outputs[1] = elcnt;
  outputs[2] = planar;

  free( Spinepath );
  free( Alefts );
  free( Arights );
  free( Elefts );
  free( Erights );
  free( alfirst );
  free( arfirst );
  free( erfirst );
  free( elfirst );
  free( arnext );

  free( alnext );
  free( ernext );
  free( elnext );
  // free(nextfreear);
  // free(nextfreeal);

  // fprintf(stdout,"Aerturning2 \n");
  free( Ae );
  free( Ee );

  // fprintf(stdout,"h: %d \n",h);
}

void Embedding_Given_Alpha( int* T, int* A, int* outputs, int* cyclic_adj, int* vidx, int e0, int t, int* alpha, int* newvadj, int* neweadj, int* newvidx, int M, int* deg, int* parent, int* backedge, int* newel, int N )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Embedding_Given_Alpha                                                  */
  /*                                                                         */
  /*  Creates an embedding from the output of the planarity check algorithm. */
  /*                                                                         */
  /***************************************************************************/

  int acnt = 0;
  int tcnt = 0;
  int x;
  int y;
  int i;
  int j;
  int r;
  int s;
  if ( backedge[e0 - 1] == 0 )
  {
    x = newel[e0 - 1];
    y = newel[M + e0 - 1];
  }
  else
  {
    x = newel[M + e0 - 1];
    y = newel[e0 - 1];
  }
  int* Spinepath;
  Spinepath = (int*)calloc( N, sizeof( int ) );
  int currv = y;
  int curre = e0;
  int scnt = 0;
  int w0;
  if ( backedge[e0 - 1] == 0 )
  {
    for ( i = 1; i <= M; i++ )
    {
      if ( backedge[curre - 1] == 1 )
      {
        w0 = currv;
        break;
      }
      scnt++;
      Spinepath[scnt - 1] = currv;
      curre = neweadj[newvidx[currv - 1]];
      currv = newvadj[newvidx[currv - 1]];
    }
  }
  else
  {
    curre = e0;
    w0 = y;
  }

  int* Tdash;
  int* Adash;
  int* AL;
  int* AR;
  Tdash = (int*)calloc( N, sizeof( int ) );
  Adash = (int*)calloc( N, sizeof( int ) );
  AL = (int*)calloc( N, sizeof( int ) );
  AR = (int*)calloc( N, sizeof( int ) );
  int tdashcnt = 0;
  int adashcnt = 0;
  int alcnt = 0;
  int arcnt = 0;
  int ta;
  int ej;
  int wj;
  int wjminus;
  int ejminus;
  int Ldash;
  int* op;
  op = (int*)calloc( N, sizeof( int ) );
  int opcnt = 0;
  int outputs2[2];

  T[0] = neweadj[newvidx[Spinepath[scnt - 1] - 1]];
  tcnt = 1;
  for ( j = scnt; j >= 1; j-- )
  {
    wj = Spinepath[j - 1];

    for ( r = 2; r <= deg[wj - 1]; r++ )
    {
      ej = neweadj[newvidx[wj - 1] + r - 1];
      if ( backedge[ej - 1] == 0 )
      {
        if ( alpha[ej - 1] == -1 )
        {
          if ( t == -1 )
          {
            ta = 1;
          }
          else
          {
            ta = -1;
          }
        }
        else
        {
          if ( t == 1 )
          {
            ta = 1;
          }
          else
          {
            ta = -1;
          }
        }
        Embedding_Given_Alpha( Tdash, Adash, outputs2, cyclic_adj, vidx, ej, ta, alpha, newvadj, neweadj, newvidx, M, deg, parent, backedge, newel, N );
        tdashcnt = outputs2[0];
        adashcnt = outputs2[1];
      }
      else
      {
        Tdash[0] = ej;
        Adash[0] = ej;
        tdashcnt = 1;
        adashcnt = 1;
      }
      if ( t == alpha[ej - 1] )
      {
        for ( s = tcnt; s >= 1; s-- )
        {
          T[s + tdashcnt - 1] = T[s - 1];
        }
        for ( s = 1; s <= tdashcnt; s++ )
        {
          T[s - 1] = Tdash[s - 1];
        }
        tcnt = tcnt + tdashcnt;
        for ( s = 1; s <= adashcnt; s++ )
        {
          AL[alcnt + s - 1] = Adash[s - 1];
        }
        alcnt = alcnt + adashcnt;
      }
      else
      {
        for ( s = 1; s <= tdashcnt; s++ )
        {
          T[tcnt + s - 1] = Tdash[s - 1];
        }
        tcnt = tcnt + tdashcnt;
        for ( s = arcnt; s >= 1; s-- )
        {
          AR[s + adashcnt - 1] = AR[s - 1];
        }
        for ( s = 1; s <= adashcnt; s++ )
        {
          AR[s - 1] = Adash[s - 1];
        }
        arcnt = arcnt + adashcnt;
      }
    }

    wjminus = parent[wj - 1];

    for ( r = newvidx[wjminus - 1] + 1; r <= newvidx[wjminus]; r++ )
    {
      if ( newvadj[r - 1] == wj )
      {
        ejminus = neweadj[r - 1];
        break;
      }
    }

    op[0] = ejminus;
    for ( r = 1; r <= tcnt; r++ )
    {
      op[r] = T[r - 1];
    }

    opcnt = 0;
    for ( r = vidx[wj - 1] + 1; r <= vidx[wj]; r++ )
    { // old undirected vidx here.
      opcnt++;
      cyclic_adj[r - 1] = op[opcnt - 1];
    }
    tcnt = 0;
    for ( r = 1; r <= alcnt; r++ )
    {
      if ( newel[AL[r - 1] - 1] == wjminus || newel[AL[r - 1] + M - 1] == wjminus )
      {
        tcnt++;
        T[tcnt - 1] = AL[r - 1];
      }
    }
    tcnt++;
    T[tcnt - 1] = ejminus;
    for ( r = 1; r <= arcnt; r++ )
    {
      if ( newel[AR[r - 1] - 1] == wjminus || newel[AR[r - 1] + M - 1] == wjminus )
      {
        tcnt++;
        T[tcnt - 1] = AR[r - 1];
      }
    }
  }
  acnt = 0;
  for ( r = 1; r <= arcnt; r++ )
  {
    acnt++;
    A[acnt - 1] = AR[r - 1];
  }
  acnt++;
  A[acnt - 1] = neweadj[newvidx[Spinepath[scnt - 1] - 1]];
  for ( r = 1; r <= alcnt; r++ )
  {
    acnt++;
    A[acnt - 1] = AL[r - 1];
  }
  outputs[0] = tcnt;
  outputs[1] = acnt;
  free( Spinepath );
  free( Tdash );
  free( Adash );
  free( AL );
  free( AR );
  free( op );
}

int Embedding_Given_Crossings( int N, int M, int* origclab, int* origcidx, int* elorig, int min_scheme, int stop_crossings, int bigface_depth, int verbose, int* finClab, int* finCidx )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Embedding_Given_Crossings                                              */
  /*                                                                         */
  /*  Constructs a valid embedding from a set of edge crossings. First, the  */
  /*  crossings are turned into vertices. If subdivisions are required this  */
  /*  can be done multiple ways. Hence, we only permit graphs which have     */
  /*  already been subdivided. Next, we perform a planarity check on the new */
  /*  graph to ensure the crossing list is valid. Then, the data from that   */
  /*  check permits us to produce an embedding.                              */
  /*                                                                         */
  /*  The planarity check algorithm is due to Tarjan and Hopcroft.           */
  /*                                                                         */
  /***************************************************************************/
  int i;
  int j;
  int k;
  int r;
  int tt;
  int ce;
  int* el;
  for ( i = 1; i <= M; i++ )
  {
    for ( j = origcidx[i - 1] + 1; j <= origcidx[i]; j++ )
    {
      for ( k = j + 1; k <= origcidx[i]; k++ )
      {
        if ( origclab[k - 1] == origclab[j - 1] )
        {
          fprintf( stdout, "multi-edge crossings. These edges need subdivided in order to solve.\n" );
          return ( -1 );
        }
      }
    }
  }
  int Mextended = M + origcidx[M];
  el = (int*)calloc( 2 * Mextended, sizeof( int ) );
  for ( i = 1; i <= M; i++ )
  {
    el[i - 1] = elorig[i - 1];
    el[Mextended + i - 1] = elorig[M + i - 1];
  }
  int origM = M;
  int origN = N;
  int* degtest = (int*)calloc( N, sizeof( int ) );
  int* vinvmap1 = (int*)calloc( N, sizeof( int ) );
  int* einvmap1 = (int*)calloc( M, sizeof( int ) );
  for ( i = 1; i <= M; i++ )
  {
    einvmap1[i - 1] = i;
  }
  for ( i = 1; i <= N; i++ )
  {
    vinvmap1[i - 1] = i;
  }
  // delete degree one verts
  int checkdel = 1;
  while ( checkdel == 1 )
  {
    checkdel = 0;
    for ( r = 0; r < N; r++ )
    {
      degtest[r] = 0;
    }
    for ( r = 1; r <= M; r++ )
    {
      degtest[el[r - 1] - 1]++;
      degtest[el[Mextended + r - 1] - 1]++;
    }
    int dele = 0;
    int delv = 0;
    for ( i = 1; i <= M; i++ )
    {
      if ( degtest[el[i - 1] - 1] == 1 )
      {
        dele = i;
        delv = el[i - 1];
        break;
      }
      if ( degtest[el[Mextended + i - 1] - 1] == 1 )
      {
        dele = i;
        delv = el[Mextended + i - 1];
        break;
      }
    }
    if ( dele > 0 )
    {
      checkdel = 0;
      for ( i = dele; i <= M - 1; i++ )
      {
        el[i - 1] = el[i];
        el[Mextended + i - 1] = el[Mextended + i];
      }
      for ( i = delv; i <= N - 1; i++ )
      {
        vinvmap1[i - 1] = vinvmap1[i];
      }
      for ( i = dele; i <= M - 1; i++ )
      {
        einvmap1[i - 1] = einvmap1[i];
      }
      M--;
      N--;
      for ( j = 1; j <= M; j++ )
      {
        if ( el[i - 1] > delv )
        {
          el[i - 1]--;
        }
        if ( el[Mextended + i - 1] > delv )
        {
          el[Mextended + i - 1]--;
        }
      }
    }
  }
  free( degtest );
  // if given clabels had degree one crossings, those need removed.
  int* emap1;
  emap1 = (int*)calloc( origM, sizeof( int ) );
  for ( i = 1; i <= M; i++ )
  {
    emap1[einvmap1[i - 1] - 1] = i;
  }

  int csize = 0;
  for ( i = 1; i <= origM; i++ )
  {
    if ( emap1[i - 1] == 0 )
    { // then this edge was deleted so delete block
      csize = origcidx[i] - origcidx[i - 1];
      if ( csize > 0 )
      {
        for ( j = origcidx[i] + 1; j <= origcidx[origM]; j++ )
        {
          origclab[j - 1 - csize] = origclab[j - 1];
        }
        for ( j = i + 1; j <= origM + 1; j++ )
        {
          origcidx[j - 1] = origcidx[j - 1] - csize;
        }
      }
    }
    for ( j = origcidx[i - 1] + 1; j <= origcidx[i]; j++ )
    {
      if ( emap1[origclab[j - 1] - 1] == 0 )
      { // then this edge was deleted so remove single entry.
        for ( k = j; k <= origcidx[origM] - 1; k++ )
        {
          origclab[k - 1] = origclab[k];
        }
        for ( k = i + 1; k <= origM + 1; k++ )
        {
          origcidx[k - 1]--;
        }
      }
    }
  }
  int* eldel1;
  eldel1 = (int*)calloc( 2 * M, sizeof( int ) );
  for ( i = 1; i <= M; i++ )
  {
    eldel1[i - 1] = el[i - 1];
    eldel1[M + i - 1] = el[Mextended + i - 1];
  }

  int* planvlabels;
  int* elistorig;
  int* eidxorig;
  int* origedgelabel;
  int* planedges;
  planvlabels = (int*)calloc( origcidx[M], sizeof( int ) );
  elistorig = (int*)calloc( M + origcidx[M], sizeof( int ) );
  eidxorig = (int*)calloc( M + 1, sizeof( int ) );
  planedges = (int*)calloc( origcidx[M], sizeof( int ) );
  int ccross = origcidx[M] / 2;
  for ( i = 0; i <= M; i++ )
  {
    eidxorig[i] = i + origcidx[i];
  }
  for ( i = 1; i <= M; i++ )
  {
    elistorig[eidxorig[i - 1]] = i;
  }
  origedgelabel = (int*)calloc( M + origcidx[M], sizeof( int ) );
  for ( i = 1; i <= M; i++ )
  {
    origedgelabel[i - 1] = i;
  }
  int Ndel1 = N;
  int Mdel1 = M;
  int f1;
  for ( i = 1; i <= Mdel1; i++ )
  {
    if ( origcidx[i] - origcidx[i - 1] != 0 )
    {
      tt = el[Mextended + i - 1];
    }
    for ( j = 1; j <= origcidx[i] - origcidx[i - 1]; j++ )
    {
      M++;
      elistorig[eidxorig[i - 1] + j] = M;
      origedgelabel[M - 1] = i;
      if ( planvlabels[origcidx[i - 1] + j - 1] == 0 )
      {
        N++;
        planvlabels[origcidx[i - 1] + j - 1] = N;
        ce = origclab[origcidx[i - 1] + j - 1];
        for ( r = origcidx[ce - 1]; r < origcidx[ce]; r++ )
        {
          if ( origclab[r] == i )
          {
            f1 = r;
            break;
          }
        }
        planvlabels[r] = N;
        planedges[N - Ndel1 - 1] = i;
        planedges[ccross + N - Ndel1 - 1] = origclab[origcidx[i - 1] + j - 1];
      }
    }
    if ( origcidx[i] - origcidx[i - 1] > 0 )
    {
      el[M - 1] = tt;
      el[Mextended + M - 1] = planvlabels[origcidx[i] - 1];
      el[Mextended + i - 1] = planvlabels[origcidx[i - 1]];
    }
    for ( j = 1; j <= origcidx[i] - origcidx[i - 1] - 1; j++ )
    {
      if ( planvlabels[origcidx[i - 1] + j - 1] < planvlabels[origcidx[i - 1] + j] )
      {
        el[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[origcidx[i - 1] + j - 1];
        el[Mextended + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[origcidx[i - 1] + j];
      }
      else
      {
        el[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[origcidx[i - 1] + j];
        el[Mextended + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[origcidx[i - 1] + j - 1];
      }
    }
  }
  int M1 = Mextended;
  Mextended = M;
  // planarized info
  int* Sixidxtotal;
  int* incidentetotal;
  Sixidxtotal = (int*)calloc( 6 * ( M + 20 ), sizeof( int ) );
  incidentetotal = (int*)calloc( N, sizeof( int ) );
  int Mlimit = M + 20;
  for ( i = 1; i <= M; i++ )
  {
    Sixidxtotal[i - 1] = el[i - 1];
    Sixidxtotal[Mlimit + i - 1] = el[M + i - 1];
  }

  int* deg;
  int* vadj;
  int* eadj;
  int* vidx;
  deg = (int*)calloc( N, sizeof( int ) );
  vadj = (int*)calloc( 2 * M, sizeof( int ) );
  eadj = (int*)calloc( 2 * M, sizeof( int ) );
  vidx = (int*)calloc( N + 1, sizeof( int ) );
  for ( i = 1; i <= M; i++ )
  {
    deg[el[i - 1] - 1]++;
    deg[el[M + i - 1] - 1]++;
  }
  for ( i = 1; i <= N; i++ )
  {
    vidx[i] = vidx[i - 1] + deg[i - 1];
  }
  int* degidx;
  degidx = (int*)calloc( N, sizeof( int ) );
  int v1;
  int v2;
  for ( i = 1; i <= M; i++ )
  {
    v1 = el[i - 1];
    v2 = el[Mdel1 + origcidx[Mdel1] + i - 1];
    degidx[v1 - 1]++;
    degidx[v2 - 1]++;
    vadj[vidx[v1 - 1] + degidx[v1 - 1] - 1] = v2;
    vadj[vidx[v2 - 1] + degidx[v2 - 1] - 1] = v1;
    eadj[vidx[v1 - 1] + degidx[v1 - 1] - 1] = i;
    eadj[vidx[v2 - 1] + degidx[v2 - 1] - 1] = i;
  }

  int maxdeg = 0;
  for ( i = 0; i < N; i++ )
  {
    if ( deg[i] > maxdeg )
    {
      maxdeg = deg[i];
    }
  }
  int Mplan = M;
  int* visited;
  int* clst;
  int* ccidx;
  int* stack;
  visited = (int*)calloc( N, sizeof( int ) );
  clst = (int*)calloc( N, sizeof( int ) );
  ccidx = (int*)calloc( N, sizeof( int ) );
  stack = (int*)calloc( N, sizeof( int ) );

  int ctcnt = 0;
  int ccnt = 0;
  int currv;

  for ( i = 1; i <= N; i++ )
  {
    if ( visited[i - 1] == 0 )
    {
      ctcnt++;
      ccnt++;
      clst[ccnt - 1] = i;
      int stackpos = 0;
      int stackmax = 1;
      stack[0] = i;
      visited[i - 1] = 1;
      while ( stackpos < stackmax )
      {
        stackpos++;
        currv = stack[stackpos - 1];
        for ( j = vidx[currv - 1] + 1; j <= vidx[currv]; j++ )
        {
          if ( visited[vadj[j - 1] - 1] == 0 )
          {
            visited[vadj[j - 1] - 1] = 1;
            ccnt++;
            clst[ccnt - 1] = vadj[j - 1];
            stackmax++;
            stack[stackmax - 1] = vadj[j - 1];
          }
        }
      }
      ccidx[ctcnt] = ccnt;
    }
  }

  int z;
  int currN;
  int vert;
  int Mcomplimit;
  int e1;
  int e2;
  int curre;
  int acompcnt;
  int tcount;
  int ecnt;
  int currM;
  // fprintf(stdout,"h1 \n");
  for ( z = 1; z <= ctcnt; z++ )
  {
    if ( ccidx[z] - ccidx[z - 1] < 5 )
    {
      continue;
    }
    currN = ccidx[z] - ccidx[z - 1];
    int* LP;
    int* vinvmapcomp;
    int* einvmapcomp;
    int* vmapcomp;
    int* visited2;
    LP = (int*)calloc( currN, sizeof( int ) );
    vinvmapcomp = (int*)calloc( currN, sizeof( int ) );
    einvmapcomp = (int*)calloc( M, sizeof( int ) );
    vmapcomp = (int*)calloc( N, sizeof( int ) );
    visited2 = (int*)calloc( currN, sizeof( int ) );
    for ( i = 0; i < currN; i++ )
    {
      vinvmapcomp[i] = clst[ccidx[z - 1] + i];
    }
    for ( i = 0; i < currN; i++ )
    {
      vmapcomp[vinvmapcomp[i] - 1] = i + 1;
    }
    int* currdeg;
    int* currvadj;
    int* curreadj;
    int* currvidx;
    currM = 0;
    // currM needs calculated first for allocations.
    for ( i = 1; i <= currN; i++ )
    {
      vert = vinvmapcomp[i - 1];
      for ( j = vidx[vert - 1] + 1; j <= vidx[vert]; j++ )
      {
        if ( vert < vadj[j - 1] )
        {
          currM++;
        }
      }
    }
    // fprintf(stdout,"h2 \n");
    int* elcomp;
    currdeg = (int*)calloc( currN, sizeof( int ) );
    currvadj = (int*)calloc( 2 * currM, sizeof( int ) );
    curreadj = (int*)calloc( 2 * currM, sizeof( int ) );
    currvidx = (int*)calloc( currN + 1, sizeof( int ) );
    elcomp = (int*)calloc( 2 * currM, sizeof( int ) );
    currM = 0;
    for ( i = 1; i <= currN; i++ )
    {
      vert = vinvmapcomp[i - 1];
      for ( j = vidx[vert - 1] + 1; j <= vidx[vert]; j++ )
      {
        if ( vert < vadj[j - 1] )
        {
          currM++;
          if ( i < vmapcomp[vadj[j - 1] - 1] )
          {
            elcomp[currM - 1] = i;
            elcomp[M + currM - 1] = vmapcomp[vadj[j - 1] - 1];
          }
          else
          {
            elcomp[M + currM - 1] = i;
            elcomp[currM - 1] = vmapcomp[vadj[j - 1] - 1];
          }
          einvmapcomp[currM - 1] = eadj[j - 1];
          currdeg[i - 1]++;
          currdeg[vmapcomp[vadj[j - 1] - 1] - 1]++;
        }
      }
    }
    // fprintf(stdout,"h3 \n");
    for ( i = 1; i <= currN; i++ )
    {
      currvidx[i] = currvidx[i - 1] + currdeg[i - 1];
    }
    int* degidx;
    degidx = (int*)calloc( currN, sizeof( int ) );
    for ( i = 1; i <= currM; i++ )
    {
      v1 = elcomp[i - 1];
      v2 = elcomp[M + i - 1];
      degidx[v1 - 1]++;
      degidx[v2 - 1]++;
      currvadj[currvidx[v1 - 1] + degidx[v1 - 1] - 1] = v2;
      currvadj[currvidx[v2 - 1] + degidx[v2 - 1] - 1] = v1;
      curreadj[currvidx[v1 - 1] + degidx[v1 - 1] - 1] = i;
      curreadj[currvidx[v2 - 1] + degidx[v2 - 1] - 1] = i;
    }

    // now need copies for extensions to biconnected component.
    // these need to be larger to fit the potential new edges.
    Mcomplimit = currM + 20;
    int* currvadj2 = (int*)calloc( 2 * Mcomplimit, sizeof( int ) );
    int* curreadj2 = (int*)calloc( 2 * Mcomplimit, sizeof( int ) );
    int* currvidx2 = (int*)calloc( currN + 1, sizeof( int ) );
    int* elcomp2 = (int*)calloc( 2 * Mcomplimit, sizeof( int ) );

    for ( i = 1; i <= 2 * currM; i++ )
    {
      currvadj2[i - 1] = currvadj[i - 1];
      curreadj2[i - 1] = curreadj[i - 1];
    }
    for ( i = 1; i <= currM; i++ )
    {
      elcomp2[i - 1] = elcomp[i - 1];
      elcomp2[Mcomplimit + i - 1] = elcomp[currM + i - 1];
    }
    for ( i = 1; i <= currN + 1; i++ )
    {
      currvidx2[i - 1] = currvidx[i - 1];
    }

    int* Tcheck = (int*)calloc( currM, sizeof( int ) );
    int* T = (int*)calloc( currM, sizeof( int ) );
    int* p = (int*)calloc( currN, sizeof( int ) );
    int* clst2 = (int*)calloc( 2 * currN, sizeof( int ) );
    int* ccidx2 = (int*)calloc( 2 * currN, sizeof( int ) );
    int* elst = (int*)calloc( currM, sizeof( int ) );
    int* eidx = (int*)calloc( currM, sizeof( int ) );
    int* artvert = (int*)calloc( currN, sizeof( int ) );
    int artcnt = 0;
    int outputs[8];

    ccnt = 0;
    tcount = 0;
    acompcnt = 0;
    ecnt = 0;
    int M2 = M;

    UDFS5( visited2, LP, 1, Tcheck, T, tcount, p, clst2, ccnt, ccidx2, acompcnt, ecnt, elst, eidx, artvert, currM, currvadj2, currvidx2, elcomp2, curreadj2, 1, elcomp, currvadj, currvidx, curreadj, outputs, Mcomplimit, artcnt, currN, currM );

    M2 = outputs[5];
    Mcomplimit = outputs[6];
    artcnt = outputs[7];

    for ( i = 1; i <= artcnt; i++ )
    {
      if ( vinvmapcomp[artvert[i - 1] - 1] > Ndel1 )
      { // then this art vert came from a crossing and can be removed.
        e1 = planedges[vinvmapcomp[artvert[i - 1] - 1] - Ndel1 - 1];
        e2 = planedges[ccross + vinvmapcomp[artvert[i - 1] - 1] - Ndel1 - 1];
        for ( j = origcidx[e1 - 1] + 1; j <= origcidx[e1]; j++ )
        {
          if ( planvlabels[j - 1] == vinvmapcomp[artvert[i - 1] - 1] )
          {
            for ( k = j; k <= origcidx[Mdel1] - 1; k++ )
            {
              origclab[k - 1] = origclab[k];
              planvlabels[k - 1] = planvlabels[k];
            }
            for ( k = e1 + 1; k <= Mdel1; k++ )
            {
              origcidx[k - 1]--;
            }
            break;
          }
        }
        for ( j = origcidx[e2 - 1] + 1; j <= origcidx[e2]; j++ )
        {
          if ( planvlabels[j - 1] == vinvmapcomp[artvert[i - 1] - 1] )
          {
            for ( k = j; k <= origcidx[Mdel1] - 1; k++ )
            {
              origclab[k - 1] = origclab[k];
              planvlabels[k - 1] = planvlabels[k];
            }
            for ( k = e2 + 1; k <= Mdel1; k++ )
            {
              origcidx[k - 1]--;
            }
            break;
          }
        }
      }
    }
    int* currdeg2;
    currdeg2 = (int*)calloc( currN, sizeof( int ) );
    for ( i = 1; i <= currN; i++ )
    {
      for ( j = currvidx2[i - 1] + 1; j <= currvidx2[i]; j++ )
      {
        currdeg2[i - 1]++;
      }
      p[i - 1] = 0;
      LP[i - 1] = 0;
    }
    // reset LP and p

    int* vmapdfs = (int*)calloc( currN, sizeof( int ) );
    int* lowpt2 = (int*)calloc( currN, sizeof( int ) );
    int* backedge = (int*)calloc( M2, sizeof( int ) );

    DFSrelabel4( vmapdfs, LP, 1, p, 0, lowpt2, backedge, curreadj2, 1, currvidx2, currvadj2, currN, outputs );

    int* newel;
    newel = (int*)calloc( 2 * M2, sizeof( int ) );
    for ( i = 1; i <= M2; i++ )
    {
      if ( vmapdfs[elcomp2[i - 1] - 1] < vmapdfs[elcomp2[Mcomplimit + i - 1] - 1] )
      {
        newel[i - 1] = vmapdfs[elcomp2[i - 1] - 1];
        newel[M2 + i - 1] = vmapdfs[elcomp2[Mcomplimit + i - 1] - 1];
      }
      else
      {
        newel[M2 + i - 1] = vmapdfs[elcomp2[i - 1] - 1];
        newel[i - 1] = vmapdfs[elcomp2[Mcomplimit + i - 1] - 1];
      }
    }
    int* buckets;
    buckets = (int*)calloc( M2, sizeof( int ) );
    int* bucketcnt;
    bucketcnt = (int*)calloc( 2 * currN + 1, sizeof( int ) );
    int* p2;
    p2 = (int*)calloc( currN, sizeof( int ) );
    int* LP2;
    LP2 = (int*)calloc( currN, sizeof( int ) );
    int* lowpt22;
    lowpt22 = (int*)calloc( currN, sizeof( int ) );
    int* invvmapdfs;
    invvmapdfs = (int*)calloc( currN, sizeof( int ) );
    for ( i = 1; i <= currN; i++ )
    {
      if ( i != 1 )
      {
        p2[vmapdfs[i - 1] - 1] = vmapdfs[p[i - 1] - 1];
      }
      LP2[vmapdfs[i - 1] - 1] = LP[i - 1];
      lowpt22[vmapdfs[i - 1] - 1] = lowpt2[i - 1];
      invvmapdfs[vmapdfs[i - 1] - 1] = i;
    }
    for ( i = 0; i < currN; i++ )
    {
      p[i] = p2[i];
      lowpt2[i] = lowpt22[i];
      LP[i] = LP2[i];
    }
    int c;
    int* bucketfirst;
    bucketfirst = (int*)calloc( 2 * currN + 1, sizeof( int ) );
    int* bucketnext;
    bucketnext = (int*)calloc( M2, sizeof( int ) );
    int freeb = 1;
    int last = 0;
    int qv;
    for ( i = 1; i <= M2; i++ )
    {
      if ( backedge[i - 1] == 1 )
      {
        c = 2 * newel[i - 1];
      }
      else
      {
        if ( p[newel[i - 1] - 1] == newel[M2 + i - 1] )
        {
          if ( lowpt2[newel[i - 1] - 1] >= newel[M2 + i - 1] )
          {
            c = 2 * LP[newel[i - 1] - 1];
          }
          else
          {
            c = 2 * LP[newel[i - 1] - 1] + 1;
          }
        }
        else
        {
          if ( lowpt2[newel[M2 + i - 1] - 1] >= newel[i - 1] )
          {
            c = 2 * LP[newel[M2 + i - 1] - 1];
          }
          else
          {
            c = 2 * LP[newel[M2 + i - 1] - 1] + 1;
          }
        }
      }

      if ( bucketfirst[c - 1] == 0 )
      {
        bucketfirst[c - 1] = freeb;
      }
      else
      {
        qv = bucketfirst[c - 1];
        for ( j = 1; j <= N; j++ )
        {
          if ( qv == 0 )
          {
            bucketnext[last - 1] = freeb;
            break;
          }
          last = qv;
          qv = bucketnext[qv - 1];
        }
      }
      buckets[freeb - 1] = i;
      freeb++;
    }

    int* newvadj;
    int* newvidx;
    int* neweadj;
    newvadj = (int*)calloc( 2 * M2, sizeof( int ) );
    neweadj = (int*)calloc( 2 * M2, sizeof( int ) );
    newvidx = (int*)calloc( currN + 1, sizeof( int ) );
    int* newdeg;
    newdeg = (int*)calloc( currN, sizeof( int ) );

    for ( i = 0; i < M2; i++ )
    {
      if ( backedge[i] == 0 )
      {
        newdeg[newel[i] - 1]++;
      }
      else
      {
        newdeg[newel[M2 + i] - 1]++;
      }
    }
    for ( i = 0; i < currN; i++ )
    {
      newvidx[i + 1] = newvidx[i] + newdeg[i];
    }
    for ( i = 0; i < currN; i++ )
    {
      degidx[i] = 0;
    }

    for ( i = 1; i <= 2 * currN; i++ )
    {
      currv = bucketfirst[i - 1];
      for ( j = 1; j <= M2; j++ )
      {
        if ( currv == 0 )
        {
          break;
        }

        e1 = newel[buckets[currv - 1] - 1];
        e2 = newel[buckets[currv - 1] - 1 + M2];

        if ( backedge[buckets[currv - 1] - 1] == 1 )
        {
          tt = e1;
          e1 = e2;
          e2 = tt;
        }
        degidx[e1 - 1]++;
        neweadj[newvidx[e1 - 1] + degidx[e1 - 1] - 1] = buckets[currv - 1];
        newvadj[newvidx[e1 - 1] + degidx[e1 - 1] - 1] = e2;
        currv = bucketnext[currv - 1];
      }
    }
    int* alpha;
    alpha = (int*)calloc( M2, sizeof( int ) );
    int* L = (int*)calloc( currM, sizeof( int ) );
    int* EL = (int*)calloc( currM, sizeof( int ) );

    Strongly_Planar( L, EL, backedge, 1, p2, neweadj, newvadj, newvidx, newel, M2, currN, newdeg, LP2, alpha, outputs );

    if ( outputs[2] == 0 )
    {
      fprintf( stdout, "Provided crossings define a non-planar graph.\n" );
      free( el );
      free( emap1 );
      free( planvlabels );
      free( elistorig );
      free( eidxorig );
      free( planedges );
      free( origedgelabel );
      free( Sixidxtotal );
      free( incidentetotal );
      free( deg );
      free( vadj );
      free( eadj );
      free( vidx );
      free( degidx );
      free( visited );
      free( clst );
      free( ccidx );
      free( stack );

      free( LP );
      free( vinvmapcomp );
      free( einvmapcomp );
      free( vmapcomp );
      free( visited2 );
      free( currdeg );
      free( currvadj );
      free( curreadj );
      free( currvidx );
      free( elcomp );
      free( currvadj2 );
      free( curreadj2 );
      free( currvidx2 );
      free( elcomp2 );
      free( Tcheck );
      free( T );
      free( p );
      free( clst2 );
      free( ccidx2 );
      free( elst );
      free( eidx );
      free( artvert );
      free( currdeg2 );
      free( vmapdfs );
      free( lowpt2 );
      free( backedge );
      free( newel );
      free( buckets );
      free( bucketcnt );
      free( p2 );
      free( LP2 );
      free( lowpt22 );
      free( invvmapdfs );
      free( bucketfirst );
      free( bucketnext );
      free( newvadj );
      free( neweadj );
      free( newvidx );
      free( newdeg );
      free( alpha );

      free( L );
      free( EL );
      return ( -1 );
    }
    int* incidente;
    int* cyclic_adj; // which will use just the currvidx2.
    int* cyc_vidx;
    incidente = (int*)calloc( currN, sizeof( int ) );
    cyclic_adj = (int*)calloc( currvidx2[currN], sizeof( int ) );
    cyc_vidx = (int*)calloc( currN + 1, sizeof( int ) );
    for ( i = 1; i <= currN; i++ )
    {
      cyc_vidx[i] = cyc_vidx[i - 1] + currdeg2[invvmapdfs[i - 1] - 1];
    }
    for ( i = 1; i <= currM; i++ )
    {
      if ( incidente[elcomp2[i - 1] - 1] == 0 )
      {
        incidente[elcomp2[i - 1] - 1] = i;
      }
      if ( incidente[elcomp2[Mcomplimit + i - 1] - 1] == 0 )
      {
        incidente[elcomp2[Mcomplimit + i - 1] - 1] = i;
      }
    }

    int* A = (int*)calloc( currM, sizeof( int ) );

    Embedding_Given_Alpha( T, A, outputs, cyclic_adj, cyc_vidx, 1, 1, alpha, newvadj, neweadj, newvidx, M2, newdeg, p2, backedge, newel, currN );

    currv = 2;
    curre = 1;
    int scnt = 0;

    for ( i = 1; i <= M2; i++ )
    {
      if ( backedge[curre - 1] == 1 )
      {
        break;
      }
      scnt++;
      curre = neweadj[newvidx[currv - 1]];
      currv = newvadj[newvidx[currv - 1]];
    }
    for ( i = 1; i <= cyc_vidx[1]; i++ )
    {
      if ( i == 1 )
      {
        cyclic_adj[i - 1] = curre;
      }
      else
      {
        cyclic_adj[i - 1] = T[i - 2];
      }
    }

    int* Sixidx;
    Sixidx = (int*)calloc( 6 * M2, sizeof( int ) );
    int ec;
    int eb;
    int ef;

    for ( i = 1; i <= M2; i++ )
    {
      Sixidx[i - 1] = elcomp2[i - 1];
      Sixidx[M2 + i - 1] = elcomp2[Mcomplimit + i - 1];
    }
    for ( j = 1; j <= currN; j++ )
    {
      for ( k = 1; k <= currdeg2[invvmapdfs[j - 1] - 1]; k++ )
      {
        ec = cyclic_adj[cyc_vidx[j - 1] + k - 1];
        if ( k == 1 )
        {
          eb = cyclic_adj[cyc_vidx[j] - 1];
        }
        else
        {
          eb = cyclic_adj[cyc_vidx[j - 1] + k - 2];
        }
        if ( k == currdeg2[invvmapdfs[j - 1] - 1] )
        {
          ef = cyclic_adj[cyc_vidx[j - 1]];
        }
        else
        {
          ef = cyclic_adj[cyc_vidx[j - 1] + k];
        }
        if ( invvmapdfs[j - 1] == elcomp2[ec - 1] )
        {
          Sixidx[2 * M2 + ec - 1] = ef;
          Sixidx[3 * M2 + ec - 1] = eb;
        }
        else
        {
          Sixidx[4 * M2 + ec - 1] = ef;
          Sixidx[5 * M2 + ec - 1] = eb;
        }
      }
    }

    if ( M2 > M )
    {
      einvmapcomp = (int*)realloc( einvmapcomp, M2 * sizeof( int ) );
    }

    for ( i = currM + 1; i <= M2; i++ )
    {
      einvmapcomp[i - 1] = Mextended + i - currM;
    }

    if ( Mextended + ( M2 - currM ) > Mlimit )
    {
      // realloc Sixidxtotal.
      Sixidxtotal = (int*)realloc( Sixidxtotal, ( Mextended + ( M2 - currM ) + 20 ) * sizeof( int ) );
      for ( k = 5; k >= 0; k-- )
      {
        for ( j = Mlimit; j >= 1; j-- )
        {
          Sixidxtotal[k * ( Mextended + ( M2 - currM ) + 20 ) + j - 1] = Sixidxtotal[k * Mlimit + j - 1];
        }
      }
      Mlimit = Mextended + ( M2 - currM ) + 20;
    }

    for ( i = Mextended + 1; i <= Mextended + ( M2 - currM ); i++ )
    {
      if ( vinvmapcomp[elcomp2[i - 1] - 1] < vinvmapcomp[elcomp2[Mcomplimit + i - 1] - 1] )
      {
        Sixidxtotal[i - 1] = vinvmapcomp[elcomp2[i - 1] - 1];
        Sixidxtotal[Mlimit + i - 1] = vinvmapcomp[elcomp2[Mcomplimit + i - 1] - 1];
      }
      else
      {
        Sixidxtotal[i - 1] = vinvmapcomp[elcomp2[Mcomplimit + i - 1] - 1];
        Sixidxtotal[Mlimit + i - 1] = vinvmapcomp[elcomp2[i - 1] - 1];
      }
    }

    Mextended = Mextended + ( M2 - currM );

    int enew;
    int estart;
    int enewmapped;

    for ( j = 1; j <= currN; j++ )
    {
      currv = vinvmapcomp[invvmapdfs[j - 1] - 1];
      enew = incidente[invvmapdfs[j - 1] - 1];
      estart = enew;
      enewmapped = einvmapcomp[enew - 1];
      incidentetotal[currv - 1] = enewmapped;
      if ( vinvmapcomp[elcomp2[enew - 1] - 1] < vinvmapcomp[elcomp2[Mcomplimit + enew - 1] - 1] )
      {
        Sixidxtotal[2 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[2 * M2 + enew - 1] - 1];
        Sixidxtotal[3 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[3 * M2 + enew - 1] - 1];
        Sixidxtotal[4 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[4 * M2 + enew - 1] - 1];
        Sixidxtotal[5 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[5 * M2 + enew - 1] - 1];
      }
      else
      {
        Sixidxtotal[2 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[4 * M2 + enew - 1] - 1];
        Sixidxtotal[3 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[5 * M2 + enew - 1] - 1];
        Sixidxtotal[4 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[2 * M2 + enew - 1] - 1];
        Sixidxtotal[5 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[3 * M2 + enew - 1] - 1];
      }
      if ( Sixidx[enew - 1] == invvmapdfs[j - 1] )
      {
        enew = Sixidx[2 * M2 + enew - 1];
      }
      else
      {
        enew = Sixidx[4 * M2 + enew - 1];
      }
      enewmapped = einvmapcomp[enew - 1];

      for ( k = 1; k <= N; k++ )
      {
        if ( enew == estart )
        {
          break;
        }
        if ( vinvmapcomp[elcomp2[enew - 1] - 1] < vinvmapcomp[elcomp2[Mcomplimit + enew - 1] - 1] )
        {
          Sixidxtotal[2 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[2 * M2 + enew - 1] - 1];
          Sixidxtotal[3 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[3 * M2 + enew - 1] - 1];
          Sixidxtotal[4 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[4 * M2 + enew - 1] - 1];
          Sixidxtotal[5 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[5 * M2 + enew - 1] - 1];
        }
        else
        {
          Sixidxtotal[2 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[4 * M2 + enew - 1] - 1];
          Sixidxtotal[3 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[5 * M2 + enew - 1] - 1];
          Sixidxtotal[4 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[2 * M2 + enew - 1] - 1];
          Sixidxtotal[5 * Mlimit + enewmapped - 1] = einvmapcomp[Sixidx[3 * M2 + enew - 1] - 1];
        }
        if ( Sixidx[enew - 1] == invvmapdfs[j - 1] )
        {
          enew = Sixidx[2 * M2 + enew - 1];
        }
        else
        {
          enew = Sixidx[4 * M2 + enew - 1];
        }
        enewmapped = einvmapcomp[enew - 1];
      }
    }

    free( LP );
    free( vinvmapcomp );
    free( einvmapcomp );
    free( vmapcomp );
    free( visited2 );
    free( currdeg );
    free( currvadj );
    free( curreadj );
    free( currvidx );
    free( elcomp );
    free( degidx );
    free( currvadj2 );
    free( curreadj2 );
    free( currvidx2 );
    free( elcomp2 );
    free( Tcheck );
    free( T );
    free( p );
    free( clst2 );
    free( ccidx2 );
    free( elst );
    free( eidx );
    free( artvert );
    free( currdeg2 );
    free( vmapdfs );
    free( lowpt2 );
    free( backedge );
    free( newel );
    free( buckets );
    free( bucketcnt );
    free( p2 );
    free( LP2 );
    free( lowpt22 );
    free( invvmapdfs );
    free( bucketfirst );
    free( bucketnext );
    free( newvadj );
    free( neweadj );
    free( newvidx );
    free( newdeg );
    free( alpha );
    free( incidente );
    free( cyclic_adj );
    free( Sixidx );

    free( L );

    free( EL );
    free( A );
  }

  for ( i = M + 1; i <= Mextended; i++ )
  {
    e1 = Sixidxtotal[2 * Mlimit + i - 1];
    e2 = Sixidxtotal[3 * Mlimit + i - 1];
    v1 = Sixidxtotal[i - 1];
    for ( j = 1; j <= N; j++ )
    {
      if ( e1 <= Mplan )
      {
        break;
      }
      if ( Sixidxtotal[e1 - 1] == v1 )
      {
        e1 = Sixidxtotal[2 * Mlimit + e1 - 1];
      }
      else
      {
        e1 = Sixidxtotal[4 * Mlimit + e1 - 1];
      }
    }
    for ( j = 1; j <= N; j++ )
    {
      if ( e2 <= Mplan )
      {
        break;
      }
      if ( Sixidxtotal[e2 - 1] == v1 )
      {
        e2 = Sixidxtotal[3 * Mlimit + e2 - 1];
      }
      else
      {
        e2 = Sixidxtotal[5 * Mlimit + e2 - 1];
      }
    }
    if ( Sixidxtotal[e1 - 1] == v1 )
    {
      Sixidxtotal[3 * Mlimit + e1 - 1] = e2;
    }
    else
    {
      Sixidxtotal[5 * Mlimit + e1 - 1] = e2;
    }
    if ( Sixidxtotal[e2 - 1] == v1 )
    {
      Sixidxtotal[2 * Mlimit + e2 - 1] = e1;
    }
    else
    {
      Sixidxtotal[4 * Mlimit + e2 - 1] = e1;
    }
    e1 = Sixidxtotal[4 * Mlimit + i - 1];
    e2 = Sixidxtotal[5 * Mlimit + i - 1];
    v1 = Sixidxtotal[Mlimit + i - 1];
    for ( j = 1; j <= N; j++ )
    {
      if ( e1 <= Mplan )
      {
        break;
      }
      if ( Sixidxtotal[e1 - 1] == v1 )
      {
        e1 = Sixidxtotal[2 * Mlimit + e1 - 1];
      }
      else
      {
        e1 = Sixidxtotal[4 * Mlimit + e1 - 1];
      }
    }
    for ( j = 1; j <= N; j++ )
    {
      if ( e2 <= Mplan )
      {
        break;
      }
      if ( Sixidxtotal[e2 - 1] == v1 )
      {
        e2 = Sixidxtotal[3 * Mlimit + e2 - 1];
      }
      else
      {
        e2 = Sixidxtotal[5 * Mlimit + e2 - 1];
      }
    }
    if ( Sixidxtotal[e1 - 1] == v1 )
    {
      Sixidxtotal[3 * Mlimit + e1 - 1] = e2;
    }
    else
    {
      Sixidxtotal[5 * Mlimit + e1 - 1] = e2;
    }
    if ( Sixidxtotal[e2 - 1] == v1 )
    {
      Sixidxtotal[2 * Mlimit + e2 - 1] = e1;
    }
    else
    {
      Sixidxtotal[4 * Mlimit + e2 - 1] = e1;
    }
  }

  // check if any crossings have been undone by the embedding.
  int nexte1;
  int nexte2;

  for ( i = Ndel1 + 1; i <= N; i++ )
  {
    curre = incidentetotal[i - 1];
    if ( Sixidxtotal[curre - 1] == i )
    {
      nexte1 = Sixidxtotal[2 * Mlimit + curre - 1];
      nexte2 = Sixidxtotal[3 * Mlimit + curre - 1];
    }
    else
    {
      nexte1 = Sixidxtotal[4 * Mlimit + curre - 1];
      nexte2 = Sixidxtotal[5 * Mlimit + curre - 1];
    }
    if ( origedgelabel[curre - 1] == origedgelabel[nexte1 - 1] )
    {
      // fprintf(stdout,"undone crossing1 %d and %d, %d %d \n",origedgelabel[curre-1],origedgelabel[nexte2-1],curre,nexte1);
      for ( j = origcidx[origedgelabel[curre - 1] - 1] + 1; j <= origcidx[origedgelabel[curre - 1]]; j++ )
      {

        if ( origclab[j - 1] == origedgelabel[nexte2 - 1] )
        {

          for ( k = j; k <= origcidx[Mdel1] - 1; k++ )
          {
            origclab[k - 1] = origclab[k];
            planvlabels[k - 1] = planvlabels[k];
          }
          for ( k = origedgelabel[curre - 1] + 1; k <= Mdel1 + 1; k++ )
          {
            origcidx[k - 1]--;
          }
          break;
        }
      }

      for ( j = origcidx[origedgelabel[nexte2 - 1] - 1] + 1; j <= origcidx[origedgelabel[nexte2 - 1]]; j++ )
      {
        if ( origclab[j - 1] == origedgelabel[curre - 1] )
        {
          for ( k = j; k <= origcidx[Mdel1] - 1; k++ )
          {
            origclab[k - 1] = origclab[k];
            planvlabels[k - 1] = planvlabels[k];
          }
          for ( k = origedgelabel[nexte2 - 1] + 1; k <= Mdel1 + 1; k++ )
          {
            origcidx[k - 1]--;
          }
          break;
        }
      }
    }
    else if ( origedgelabel[curre - 1] == origedgelabel[nexte2 - 1] )
    {
      // fprintf(stdout,"undone crossing2 %d and %d, %d %d  \n",origedgelabel[curre-1],origedgelabel[nexte1-1],curre, nexte1);
      for ( j = origcidx[origedgelabel[curre - 1] - 1] + 1; j <= origcidx[origedgelabel[curre - 1]]; j++ )
      {
        if ( origclab[j - 1] == origedgelabel[nexte1 - 1] )
        {
          for ( k = j; k <= origcidx[Mdel1] - 1; k++ )
          {
            origclab[k - 1] = origclab[k];
            planvlabels[k - 1] = planvlabels[k];
          }
          for ( k = origedgelabel[curre - 1] + 1; k <= Mdel1 + 1; k++ )
          {
            origcidx[k - 1]--;
          }
          break;
        }
      }
      for ( j = origcidx[origedgelabel[nexte1 - 1] - 1] + 1; j <= origcidx[origedgelabel[nexte1 - 1]]; j++ )
      {
        if ( origclab[j - 1] == origedgelabel[curre - 1] )
        {
          for ( k = j; k <= origcidx[Mdel1] - 1; k++ )
          {
            origclab[k - 1] = origclab[k];
            planvlabels[k - 1] = planvlabels[k];
          }
          for ( k = origedgelabel[nexte1 - 1] + 1; k <= Mdel1 + 1; k++ )
          {
            origcidx[k - 1]--;
          }
          break;
        }
      }
    }
  }

  int first;
  int* origdlab = (int*)calloc( origcidx[Mdel1], sizeof( int ) );
  int oldv;
  int starte;
  int oe;
  int oefirst;
  int qq;

  for ( i = 1; i <= Mdel1; i++ )
  {
    for ( r = origcidx[i - 1] + 1; r <= origcidx[i]; r++ )
    {
      if ( i < origclab[r - 1] )
      {
        first = 0;
        vert = planvlabels[r - 1];
        oldv = vert;
        curre = incidentetotal[vert - 1];
        starte = curre;
        for ( j = 1; j <= 2; j++ )
        {
          if ( j == 2 )
          {
            if ( Sixidxtotal[starte - 1] == vert )
            {
              starte = Sixidxtotal[2 * Mlimit + starte - 1];
            }
            else
            {
              starte = Sixidxtotal[4 * Mlimit + starte - 1];
            }
            curre = starte;
            oldv = vert;
          }
          for ( k = 1; k <= M; k++ )
          {
            if ( Sixidxtotal[curre - 1] == oldv )
            {
              oldv = Sixidxtotal[Mlimit + curre - 1];
            }
            else
            {
              oldv = Sixidxtotal[curre - 1];
            }
            if ( oldv <= Ndel1 )
            {

              oe = origedgelabel[curre - 1];
              if ( eldel1[oe - 1] == oldv )
              {
                if ( first == 1 )
                {
                  for ( qq = origcidx[oefirst - 1] + 1; qq <= origcidx[oefirst]; qq++ )
                  {
                    if ( origclab[qq - 1] == oe )
                    {
                      origdlab[qq - 1] = -1;
                      break;
                    }
                  }
                  for ( qq = origcidx[oe - 1] + 1; qq <= origcidx[oe]; qq++ )
                  {
                    if ( origclab[qq - 1] == oefirst )
                    {
                      origdlab[qq - 1] = 1;
                      break;
                    }
                  }
                  break;
                }
                else if ( first == 2 )
                {
                  for ( qq = origcidx[oefirst - 1] + 1; qq <= origcidx[oefirst]; qq++ )
                  {
                    if ( origclab[qq - 1] == oe )
                    {
                      origdlab[qq - 1] = 1;
                      break;
                    }
                  }
                  for ( qq = origcidx[oe - 1] + 1; qq <= origcidx[oe]; qq++ )
                  {
                    if ( origclab[qq - 1] == oefirst )
                    {
                      origdlab[qq - 1] = -1;
                      break;
                    }
                  }
                  break;
                }
                oefirst = oe;
                first = 1;
                break;
              }
              else
              {
                if ( first == 1 )
                {
                  for ( qq = origcidx[oefirst - 1] + 1; qq <= origcidx[oefirst]; qq++ )
                  {
                    if ( origclab[qq - 1] == oe )
                    {
                      origdlab[qq - 1] = 1;
                      break;
                    }
                  }
                  for ( qq = origcidx[oe - 1] + 1; qq <= origcidx[oe]; qq++ )
                  {
                    if ( origclab[qq - 1] == oefirst )
                    {
                      origdlab[qq - 1] = -1;
                      break;
                    }
                  }
                  break;
                }
                else if ( first == 2 )
                {
                  for ( qq = origcidx[oefirst - 1] + 1; qq <= origcidx[oefirst]; qq++ )
                  {
                    if ( origclab[qq - 1] == oe )
                    {
                      origdlab[qq - 1] = -1;
                      break;
                    }
                  }
                  for ( qq = origcidx[oe - 1] + 1; qq <= origcidx[oe]; qq++ )
                  {
                    if ( origclab[qq - 1] == oefirst )
                    {
                      origdlab[qq - 1] = 1;
                      break;
                    }
                  }
                  break;
                }
                oefirst = oe;
                first = 2;
                break;
              }
            }
            if ( Sixidxtotal[curre - 1] == oldv )
            {
              curre = Sixidxtotal[2 * Mlimit + curre - 1];
            }
            else
            {
              curre = Sixidxtotal[4 * Mlimit + curre - 1];
            }
            if ( origedgelabel[curre - 1] != origedgelabel[starte - 1] )
            {
              for ( qq = 1; qq <= N; qq++ )
              {
                if ( origedgelabel[curre - 1] == origedgelabel[starte - 1] )
                {
                  break;
                }
                if ( Sixidxtotal[curre - 1] == oldv )
                {
                  curre = Sixidxtotal[2 * Mlimit + curre - 1];
                }
                else
                {
                  curre = Sixidxtotal[4 * Mlimit + curre - 1];
                }
              }
            }
          }
        }
      }
    }
  }

  // now info for original unplanarized graph.
  int* dego = (int*)calloc( Ndel1, sizeof( int ) );
  int* vadjo = (int*)calloc( 2 * Mdel1, sizeof( int ) );
  int* eadjo = (int*)calloc( 2 * Mdel1, sizeof( int ) );
  int* vidxo = (int*)calloc( Ndel1 + 1, sizeof( int ) );
  for ( i = 1; i <= Mdel1; i++ )
  {
    dego[eldel1[i - 1] - 1]++;
    dego[eldel1[Mdel1 + i - 1] - 1]++;
  }
  for ( i = 1; i <= Ndel1; i++ )
  {
    vidxo[i] = vidxo[i - 1] + dego[i - 1];
  }
  int* degidxo = (int*)calloc( Ndel1, sizeof( int ) );
  for ( i = 1; i <= Mdel1; i++ )
  {
    v1 = eldel1[i - 1];
    v2 = eldel1[Mdel1 + i - 1];
    degidxo[v1 - 1]++;
    degidxo[v2 - 1]++;
    vadjo[vidxo[v1 - 1] + degidxo[v1 - 1] - 1] = v2;
    vadjo[vidxo[v2 - 1] + degidxo[v2 - 1] - 1] = v1;
    eadjo[vidxo[v1 - 1] + degidxo[v1 - 1] - 1] = i;
    eadjo[vidxo[v2 - 1] + degidxo[v2 - 1] - 1] = i;
  }
  int compN;
  int startv;
  int Final_crossings = 0;
  int fcnt = 0;
  int* finClab1 = (int*)calloc( origcidx[Mdel1], sizeof( int ) );
  int* finCidx1 = (int*)calloc( origcidx[Mdel1], sizeof( int ) );
  int* finCcnt = (int*)calloc( Mdel1, sizeof( int ) );

  for ( z = 1; z <= ctcnt; z++ )
  {
    compN = 0;
    startv = 0;
    for ( i = ccidx[z - 1] + 1; i <= ccidx[z]; i++ )
    {
      if ( clst[i - 1] <= Ndel1 )
      {
        compN++;
        if ( startv == 0 )
        {
          startv = clst[i - 1];
        }
      }
    }
    int* visitedo = (int*)calloc( Ndel1, sizeof( int ) );
    int* LPo = (int*)calloc( Ndel1, sizeof( int ) );
    int* Tchecko = (int*)calloc( Mdel1, sizeof( int ) );
    int* To = (int*)calloc( Mdel1, sizeof( int ) );
    int* po = (int*)calloc( Ndel1, sizeof( int ) );
    int* clst2o = (int*)calloc( 2 * Mdel1, sizeof( int ) );
    int* ccidx2o = (int*)calloc( Ndel1, sizeof( int ) );
    int* elsto = (int*)calloc( Mdel1, sizeof( int ) );
    int* eidxo = (int*)calloc( Mdel1, sizeof( int ) );
    int* artvertso = (int*)calloc( Ndel1, sizeof( int ) );

    ccnt = 0;
    acompcnt = 0;
    tcount = 0;
    ecnt = 0;
    int outputs[3];

    UDFS4( startv, 1, visitedo, LPo, Tchecko, To, 0, po, eldel1, vidxo, vadjo, clst2o, ccnt, eadjo, Ndel1, ccidx2o, acompcnt, ecnt, elsto, eidxo, artvertso, 0, outputs, Mdel1 );
    acompcnt = outputs[3];

    for ( i = 1; i <= acompcnt; i++ )
    {
      currN = ccidx2o[i] - ccidx2o[i - 1];
      if ( currN < 5 )
      {
        continue;
      }
      currM = eidxo[i] - eidxo[i - 1];

      int* currSixidx = (int*)calloc( 6 * currM, sizeof( int ) );
      int cnt = 0;
      int* vmapfinal = (int*)calloc( Ndel1, sizeof( int ) );
      for ( j = 1; j <= currN; j++ )
      {
        vmapfinal[clst2o[ccidx2o[i - 1] + j - 1] - 1] = j;
      }
      int* emapfin = (int*)calloc( Mdel1, sizeof( int ) );
      int* einvmapfin = (int*)calloc( currM, sizeof( int ) );
      int cnttmp = 0;
      for ( j = eidxo[i - 1] + 1; j <= eidxo[i]; j++ )
      {
        cnttmp++;
        einvmapfin[cnttmp - 1] = elsto[j - 1];
      }
      for ( j = 1; j <= currM; j++ )
      {
        emapfin[elsto[eidxo[i - 1] + j - 1] - 1] = j;
      }
      int* vinvmapfin = (int*)calloc( currN, sizeof( int ) );
      cnttmp = 0;
      for ( j = ccidx2o[i - 1] + 1; j <= ccidx2o[i]; j++ )
      {
        cnttmp++;
        vinvmapfin[cnttmp - 1] = clst2o[j - 1];
      }
      int efu;
      int swap;

      for ( j = eidxo[i - 1] + 1; j <= eidxo[i]; j++ )
      {
        cnt++;
        efu = elsto[j - 1];
        // fprintf(stdout,"efu: %d\n",efu);
        if ( vmapfinal[eldel1[efu - 1] - 1] < vmapfinal[eldel1[Mdel1 + efu - 1] - 1] )
        {
          currSixidx[cnt - 1] = vmapfinal[eldel1[efu - 1] - 1];
          currSixidx[currM + cnt - 1] = vmapfinal[eldel1[Mdel1 + efu - 1] - 1];
          swap = 0;
        }
        else
        {
          currSixidx[cnt - 1] = vmapfinal[eldel1[Mdel1 + efu - 1] - 1];
          currSixidx[currM + cnt - 1] = vmapfinal[eldel1[efu - 1] - 1];
          swap = 1;
        }
        int ef = elistorig[eidxorig[efu - 1]];
        // fprintf(stdout,"ef: %d\n",ef);
        if ( Sixidxtotal[ef - 1] <= Ndel1 )
        {
          curre = Sixidxtotal[2 * Mlimit + ef - 1];
          for ( k = 1; k <= N; k++ )
          {
            if ( emapfin[origedgelabel[curre - 1] - 1] > 0 )
            {
              break;
            }
            if ( Sixidxtotal[curre - 1] == Sixidxtotal[ef - 1] )
            {
              curre = Sixidxtotal[2 * Mlimit + curre - 1];
            }
            else
            {
              curre = Sixidxtotal[4 * Mlimit + curre - 1];
            }
          }
          if ( swap == 0 )
          {
            currSixidx[2 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
          else
          {
            currSixidx[4 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
          curre = Sixidxtotal[3 * Mlimit + ef - 1];
          for ( k = 1; k <= N; k++ )
          {
            if ( emapfin[origedgelabel[curre - 1] - 1] > 0 )
            {
              break;
            }
            if ( Sixidxtotal[curre - 1] == Sixidxtotal[ef - 1] )
            {
              curre = Sixidxtotal[3 * Mlimit + curre - 1];
            }
            else
            {
              curre = Sixidxtotal[5 * Mlimit + curre - 1];
            }
          }
          if ( swap == 0 )
          {
            currSixidx[3 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
          else
          {
            currSixidx[5 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
        }

        int ef2 = elistorig[eidxorig[efu] - 1];
        if ( ef2 != ef )
        {
          if ( Sixidxtotal[ef2 - 1] <= Ndel1 )
          {
            curre = Sixidxtotal[2 * Mlimit + ef2 - 1];
            for ( k = 1; k <= N; k++ )
            {
              if ( emapfin[origedgelabel[curre - 1] - 1] > 0 )
              {
                break;
              }
              if ( Sixidxtotal[curre - 1] == Sixidxtotal[ef2 - 1] )
              {
                curre = Sixidxtotal[2 * Mlimit + curre - 1];
              }
              else
              {
                curre = Sixidxtotal[4 * Mlimit + curre - 1];
              }
            }
            if ( swap == 0 )
            {
              currSixidx[4 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
            }
            else
            {
              currSixidx[2 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
            }
            curre = Sixidxtotal[3 * Mlimit + ef2 - 1];
            for ( k = 1; k <= N; k++ )
            {
              if ( emapfin[origedgelabel[curre - 1] - 1] > 0 )
              {
                break;
              }
              if ( Sixidxtotal[curre - 1] == Sixidxtotal[ef2 - 1] )
              {
                curre = Sixidxtotal[3 * Mlimit + curre - 1];
              }
              else
              {
                curre = Sixidxtotal[5 * Mlimit + curre - 1];
              }
            }
            if ( swap == 0 )
            {
              currSixidx[5 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
            }
            else
            {
              currSixidx[3 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
            }
          }
        }
        else
        {
          curre = Sixidxtotal[4 * Mlimit + ef2 - 1];
          for ( k = 1; k <= N; k++ )
          {
            if ( emapfin[origedgelabel[curre - 1] - 1] > 0 )
            {
              break;
            }
            if ( Sixidxtotal[curre - 1] == Sixidxtotal[Mlimit + ef2 - 1] )
            {
              curre = Sixidxtotal[2 * Mlimit + curre - 1];
            }
            else
            {
              curre = Sixidxtotal[4 * Mlimit + curre - 1];
            }
          }
          if ( swap == 0 )
          {
            currSixidx[4 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
          else
          {
            currSixidx[2 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
          curre = Sixidxtotal[5 * Mlimit + ef2 - 1];
          for ( k = 1; k <= N; k++ )
          {
            if ( emapfin[origedgelabel[curre - 1] - 1] > 0 )
            {
              break;
            }
            if ( Sixidxtotal[curre - 1] == Sixidxtotal[Mlimit + ef2 - 1] )
            {
              curre = Sixidxtotal[3 * Mlimit + curre - 1];
            }
            else
            {
              curre = Sixidxtotal[5 * Mlimit + curre - 1];
            }
          }
          if ( swap == 0 )
          {
            currSixidx[5 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
          else
          {
            currSixidx[3 * currM + cnt - 1] = emapfin[origedgelabel[curre - 1] - 1];
          }
        }
      }

      int* currclab = (int*)calloc( origcidx[Mdel1], sizeof( int ) );
      int* currdlab = (int*)calloc( origcidx[Mdel1], sizeof( int ) );
      int* currcidx = (int*)calloc( currM + 1, sizeof( int ) );
      // int *currplanvlabels=(int*)calloc(origcidx[Mdel1],sizeof(int));
      int maxcrossings = 0;
      int currx;
      cnt = 0;
      for ( j = 1; j <= currM; j++ )
      {
        oe = elsto[eidxo[i - 1] + j - 1];
        // fprintf(stdout,"oe: %d \n",oe);
        currx = 0;
        if ( eldel1[oe - 1] == clst2o[ccidx2o[i - 1] + currSixidx[j - 1] - 1] )
        {
          for ( k = origcidx[oe - 1] + 1; k <= origcidx[oe]; k++ )
          {
            // fprintf(stdout,"k1: %d \n",k);
            if ( emapfin[origclab[k - 1] - 1] > 0 )
            {
              currx++;
              cnt++;
              currclab[cnt - 1] = emapfin[origclab[k - 1] - 1];
              // currplanvlabels[cnt-1]=planvlabels[k-1];
              int pe = origclab[k - 1];
              if ( eldel1[pe - 1] == vinvmapfin[currSixidx[emapfin[pe - 1] - 1] - 1] )
              {
                currdlab[cnt - 1] = origdlab[k - 1];
              }
              else
              {
                currdlab[cnt - 1] = -origdlab[k - 1];
              }
            }
          }
        }
        else
        {
          for ( k = origcidx[oe]; k >= origcidx[oe - 1] + 1; k-- )
          {
            // fprintf(stdout,"k2: %d \n",k);
            if ( emapfin[origclab[k - 1] - 1] > 0 )
            {
              currx++;
              cnt++;
              currclab[cnt - 1] = emapfin[origclab[k - 1] - 1];
              // currplanvlabels[cnt-1]=planvlabels[k-1];
              int pe = origclab[k - 1];
              if ( eldel1[pe - 1] == vinvmapfin[currSixidx[emapfin[pe - 1] - 1] - 1] )
              {
                currdlab[cnt - 1] = -origdlab[k - 1];
              }
              else
              {
                currdlab[cnt - 1] = origdlab[k - 1];
              }
            }
          }
        }
        if ( currx > maxcrossings )
        {
          maxcrossings = currx;
        }
        currcidx[j] = cnt;
      }

      int* currincidente = (int*)calloc( currN, sizeof( int ) );
      for ( j = 1; j <= currM; j++ )
      {
        if ( currincidente[currSixidx[j - 1] - 1] == 0 )
        {
          currincidente[currSixidx[j - 1] - 1] = j;
        }
        if ( currincidente[currSixidx[currM + j - 1] - 1] == 0 )
        {
          currincidente[currSixidx[currM + j - 1] - 1] = j;
        }
      }
      int current_crossings = currcidx[currM] / 2;
      int* currdego = (int*)calloc( currN, sizeof( int ) );
      for ( j = 1; j <= currM; j++ )
      {
        currdego[currSixidx[j - 1] - 1]++;
        currdego[currSixidx[currM + j - 1] - 1]++;
      }

      Quick_Cross_Main_Loop( &currincidente, &currSixidx, &currclab, &currdlab, &currcidx, currN, currM, &currdego, current_crossings, min_scheme, bigface_depth, stop_crossings, true, currN, outputs, verbose );

      currM = outputs[1];
      Final_crossings = Final_crossings + outputs[0];
      // fprintf(stdout,"subcr: %d fin cr: %d \n",outputs[0],Final_crossings);
      // currclab has now changed inside main loop.
      if ( fcnt + currcidx[currM] > origM )
      {
        finClab1 = (int*)realloc( finClab1, ( fcnt + currcidx[M] ) * sizeof( int ) );
        finCidx1 = (int*)realloc( finCidx1, ( fcnt + currcidx[M] ) * sizeof( int ) );
      }
      for ( j = 1; j <= currM; j++ )
      {
        if ( eldel1[einvmapfin[j - 1] - 1] == clst2o[ccidx2o[i - 1] + currSixidx[j - 1] - 1] )
        {
          for ( k = currcidx[j - 1] + 1; k <= currcidx[j]; k++ )
          {
            fcnt++;
            finClab1[fcnt - 1] = einvmap1[einvmapfin[currclab[k - 1] - 1] - 1];
            finCidx1[fcnt - 1] = einvmap1[einvmapfin[j - 1] - 1];
            finCcnt[einvmap1[einvmapfin[j - 1] - 1] - 1]++;
          }
        }
        else
        {
          for ( k = currcidx[j]; k >= currcidx[j - 1] + 1; k-- )
          {
            fcnt++;
            finClab1[fcnt - 1] = einvmap1[einvmapfin[currclab[k - 1] - 1] - 1];
            finCidx1[fcnt - 1] = einvmap1[einvmapfin[j - 1] - 1];
            finCcnt[einvmap1[einvmapfin[j - 1] - 1] - 1]++;
          }
        }
      }
      free( currSixidx );
      free( vmapfinal );
      free( emapfin );
      free( einvmapfin );
      free( vinvmapfin );
      free( currclab );
      free( currdlab );
      free( currcidx );
      free( currincidente );
      free( currdego );
    }
    free( visitedo );
    free( LPo );
    free( Tchecko );
    free( To );
    free( po );
    free( clst2o );
    free( ccidx2o );
    free( elsto );
    free( eidxo );
    free( artvertso );
  }
  int* fdn;
  fdn = (int*)calloc( origM, sizeof( int ) );
  for ( i = 1; i <= origM; i++ )
  {
    finCidx[i] = finCidx[i - 1] + finCcnt[i - 1];
  }
  for ( i = 1; i <= fcnt; i++ )
  {
    fdn[finCidx1[i - 1] - 1]++;
    finClab[finCidx[finCidx1[i - 1] - 1] - 1 + fdn[finCidx1[i - 1] - 1]] = finClab1[i - 1];
  }

  free( el );
  free( emap1 );
  free( planvlabels );
  free( elistorig );
  free( eidxorig );
  free( planedges );
  free( origedgelabel );
  free( Sixidxtotal );
  free( incidentetotal );
  free( deg );
  free( vadj );
  free( eadj );
  free( vidx );
  free( degidx );
  free( visited );
  free( clst );
  free( ccidx );
  free( stack );

  free( origdlab );
  free( dego );
  free( vadjo );
  free( eadjo );
  free( vidxo );
  free( degidxo );

  free( vinvmap1 );
  free( einvmap1 );
  free( finCcnt );

  free( finCidx1 );
  free( finClab1 );
  free( fdn );
  return ( Final_crossings );
}

void Biconnected_Runner( int* el, int N, int M, int min_scheme, int embed_scheme, int bigface_depth, int seed, int verbose, int output_scheme, int stop_crossings, char* output_filename, int* finoutputs, int** finClab, int* finCidx, long double* x, long double* y )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Biconnected_Runner                                                     */
  /*                                                                         */
  /*  This subroutine preprocesses the graph. It performs two functions:     */
  /*                                                                         */
  /* 1) It separates the graph into maximally 2-connected components. This   */
  /*    is vital as the algorithm expects the graph to be connected at all   */
  /*    times, including when a vertex (considered for moving) is removed.   */
  /*    For disconnected components this is simply a case of considering     */
  /*    them individually. For 1-connected components articulation vertices  */
  /*    are identified and the 2-connected components are built accordingly. */
  /*                                                                         */
  /* 2) Degree 2 vertices are removed ("smoothed") whenever possible (that   */
  /*    is, when their removal would not create a multi-edge).               */
  /*                                                                         */
  /*                                                                         */
  /*  An embedding is then found for each maximally 2-connected component,   */
  /*  and the components are submitted to Quick_Cross_Main_Loop.             */
  /*                                                                         */
  /***************************************************************************/

  bool acyclic = true;
  int* deg;
  deg = (int*)calloc( N, sizeof( int ) );
  int* vadj;
  vadj = (int*)calloc( 2 * M, sizeof( int ) );
  int* eadj;
  eadj = (int*)calloc( 2 * M, sizeof( int ) );
  int vidx[N + 1];
  vidx[0] = 0;
  int i;
  int current_crossings;
  int outputs2[3];
  for ( i = 0; i < M; i++ )
  {
    deg[el[i] - 1] = deg[el[i] - 1] + 1;
    deg[el[M + i] - 1] = deg[el[M + i] - 1] + 1;
  }
  for ( i = 0; i < N; i++ )
  {
    vidx[i + 1] = vidx[i] + deg[i];
  }
  int* degchk;
  degchk = (int*)calloc( N, sizeof( int ) );
  for ( i = 0; i < M; i++ )
  {
    vadj[vidx[el[i] - 1] + degchk[el[i] - 1]] = el[M + i];
    vadj[vidx[el[M + i] - 1] + degchk[el[M + i] - 1]] = el[i];
    eadj[vidx[el[i] - 1] + degchk[el[i] - 1]] = i + 1;
    eadj[vidx[el[M + i] - 1] + degchk[el[M + i] - 1]] = i + 1;
    degchk[el[i] - 1] = degchk[el[i] - 1] + 1;
    degchk[el[M + i] - 1] = degchk[el[M + i] - 1] + 1;
  }
  int checkdel = 1;
  int changed = 0;
  int v1;
  int j;
  int f1;
  int Mcnt = 0;
  int* einvmapdel1;
  einvmapdel1 = (int*)malloc( M * sizeof( int ) );
  for ( i = 0; i < M; i++ )
  {
    einvmapdel1[i] = i + 1;
  }
  int* vinvmapdel1;
  vinvmapdel1 = (int*)calloc( N, sizeof( int ) );
  for ( i = 0; i < N; i++ )
  {
    vinvmapdel1[i] = i + 1;
  }
  int* finClab1;
  int* finCcnt;
  int* finCidx1;
  finClab1 = (int*)malloc( M * sizeof( int ) );
  finCcnt = (int*)calloc( M, sizeof( int ) );
  finCidx1 = (int*)malloc( M * sizeof( int ) );
  int origM = M;
  int origN = N;
  int fcnt = 0;
  while ( checkdel == 1 )
  {
    checkdel = 0;
    for ( i = 0; i < N; i++ )
    {
      if ( deg[i] == 1 )
      {
        changed = 1;
        checkdel = 1;
        v1 = vadj[vidx[i]];
        for ( j = vidx[v1 - 1]; j < vidx[v1]; j++ )
        {
          if ( vadj[j] == i + 1 )
          {
            f1 = j;
            break;
          }
        }
        for ( j = 0; j < vidx[N]; j++ )
        {
          if ( vadj[j] > i + 1 )
          {
            vadj[j] = vadj[j] - 1;
          }
          if ( eadj[j] > eadj[vidx[i]] )
          {
            eadj[j] = eadj[j] - 1;
          }
        }
        for ( j = eadj[vidx[i]]; j < origM; j++ )
        {
          einvmapdel1[j - 1] = einvmapdel1[j];
        }
        for ( j = i + 1; j < N; j++ )
        {
          vinvmapdel1[j - 1] = vinvmapdel1[j];
        }
        for ( j = f1; j < vidx[N] - 1; j++ )
        {
          vadj[j] = vadj[j + 1];
          eadj[j] = eadj[j + 1];
        }
        deg[v1 - 1] = deg[v1 - 1] - 1;
        for ( j = v1; j <= N; j++ )
        {
          vidx[j] = vidx[j] - 1;
        }
        for ( j = vidx[i]; j < vidx[N]; j++ )
        {
          vadj[j] = vadj[j + 1];
          eadj[j] = eadj[j + 1];
        }

        for ( j = i + 1; j < N; j++ )
        {
          vidx[j] = vidx[j + 1] - 1;
        }
        for ( j = i; j < N - 1; j++ )
        {
          deg[j] = deg[j + 1];
        }
        N = N - 1;
        M = M - 1;
        break;
      }
    }
  }
  int ccnt = 0;
  int ctcnt = 0;
  int* ccidx;
  ccidx = (int*)calloc( N + 1, sizeof( int ) );
  int* visited;
  visited = (int*)calloc( N, sizeof( int ) );
  int* stack;
  stack = (int*)calloc( N, sizeof( int ) );
  int* clst;
  clst = (int*)calloc( N, sizeof( int ) );
  int stackpos;
  int stackmax;
  int currv;
  int ci;
  for ( i = 0; i < N; i++ )
  {
    if ( visited[i] == 0 )
    {
      ctcnt = ctcnt + 1;
      ccnt = ccnt + 1;
      clst[ccnt - 1] = i + 1;
      stackpos = 0;
      stackmax = 1;
      stack[0] = i + 1;
      visited[i] = 1;
      while ( stackpos < stackmax )
      {
        stackpos = stackpos + 1;
        currv = stack[stackpos - 1];
        for ( j = vidx[currv - 1]; j < vidx[currv]; j++ )
        {
          if ( visited[vadj[j] - 1] == 0 )
          {
            visited[vadj[j] - 1] = 1;
            ccnt = ccnt + 1;
            clst[ccnt - 1] = vadj[j];
            stackmax = stackmax + 1;
            stack[stackmax - 1] = vadj[j];
          }
        }
      }
      ccidx[ctcnt] = ccnt;
    }
  }

  int Final_crossings = 0;
  int Final_crossings2 = 0;
  int z;
  int k;
  for ( z = 1; z <= ctcnt; z++ )
  {

    int compN = ccidx[z] - ccidx[z - 1];
    if ( compN < 5 )
    {
      continue;
    }
    int* visited;
    visited = (int*)calloc( N, sizeof( int ) );
    int* LP;
    LP = (int*)calloc( N, sizeof( int ) );
    int* Tcheck;
    Tcheck = (int*)calloc( M, sizeof( int ) );
    int* T;
    T = (int*)calloc( M, sizeof( int ) );
    int tcount = 0;
    int* p;
    p = (int*)calloc( N, sizeof( int ) );
    int* clst2;
    clst2 = (int*)calloc( 2 * compN, sizeof( int ) );
    int* ccidx2;
    ccidx2 = (int*)calloc( compN, sizeof( int ) );
    ccnt = 0;
    int acompcnt = 0;
    int ecnt = 0;
    int* elst;
    elst = (int*)calloc( M, sizeof( int ) );
    int* eidx;
    eidx = (int*)calloc( compN, sizeof( int ) );
    int outputsA[4];

    UDFS( clst[ccidx[z - 1]], 1, visited, LP, Tcheck, T, tcount, p, el, vidx, vadj, clst2, ccnt, eadj, origN, ccidx2, acompcnt, einvmapdel1, ecnt, elst, eidx, outputsA, origM );
    acompcnt = outputsA[1];

    if ( acompcnt > 1 )
    {
      for ( i = 1; i <= acompcnt; i++ )
      {
        int relcnt = ccidx2[i] - ccidx2[i - 1];
        if ( relcnt < 5 )
        {
          continue; // no contributions
        }
        int Mcnt = eidx[i] - eidx[i - 1];
        int* currvadj;
        currvadj = (int*)calloc( 2 * Mcnt, sizeof( int ) );
        int* curreadj;
        curreadj = (int*)calloc( 2 * Mcnt, sizeof( int ) );
        int* currvidx;
        currvidx = (int*)calloc( ( relcnt + 1 ), sizeof( int ) );
        int* currdeg;
        currdeg = (int*)calloc( relcnt, sizeof( int ) );
        int* einvmaptwo;
        einvmaptwo = (int*)calloc( Mcnt, sizeof( int ) );
        for ( j = 0; j < Mcnt; j++ )
        {
          einvmaptwo[j] = elst[eidx[i - 1] + j];
        }
        int* vmaptwo;
        vmaptwo = (int*)calloc( origN, sizeof( int ) );

        int* Clabel;
        int* Dlabel;
        int* Cidx;
        int maxdeg = 0;
        Clabel = (int*)calloc( relcnt, sizeof( int ) );
        Dlabel = (int*)calloc( relcnt, sizeof( int ) );
        Cidx = (int*)calloc( ( Mcnt + 1 ), sizeof( int ) );
        int* vinvmaptwo;
        vinvmaptwo = (int*)calloc( relcnt, sizeof( int ) );
        for ( j = 0; j < relcnt; j++ )
        {
          vinvmaptwo[j] = clst2[ccidx2[i - 1] + j];
        }
        for ( j = 0; j < relcnt; j++ )
        {
          vmaptwo[vinvmaptwo[j] - 1] = j + 1;
        }
        for ( j = eidx[i - 1]; j < eidx[i]; j++ )
        {
          currdeg[vmaptwo[el[elst[j] - 1] - 1] - 1]++;
          currdeg[vmaptwo[el[origM + elst[j] - 1] - 1] - 1]++;
        }
        for ( j = 0; j < relcnt; j++ )
        {
          currvidx[j + 1] = currvidx[j] + currdeg[j];
          degchk[j] = 0;
        }
        int v1;
        int v2;
        int jc = 0;
        for ( j = eidx[i - 1]; j < eidx[i]; j++ )
        {
          jc++;
          v1 = vmaptwo[el[elst[j] - 1] - 1];
          v2 = vmaptwo[el[origM + elst[j] - 1] - 1];
          degchk[v1 - 1]++;
          degchk[v2 - 1]++;
          currvadj[currvidx[v1 - 1] + degchk[v1 - 1] - 1] = v2;
          currvadj[currvidx[v2 - 1] + degchk[v2 - 1] - 1] = v1;
          curreadj[currvidx[v1 - 1] + degchk[v1 - 1] - 1] = jc;
          curreadj[currvidx[v2 - 1] + degchk[v2 - 1] - 1] = jc;
        }

        int e1;
        int t = 1;
        int f2cnt;
        int v2chk;
        int* f2;
        f2 = (int*)malloc( relcnt * sizeof( int ) );
        int e2;
        while ( t == 1 )
        {
          t = 0;
          if ( relcnt < 5 )
          {
            break;
          }
          for ( j = 0; j < relcnt; j++ )
          {
            if ( currdeg[j] == 2 )
            {
              t = 1;
              v1 = currvadj[currvidx[j]];
              v2 = currvadj[currvidx[j] + 1];
              if ( currdeg[v2 - 1] > 2 )
              {
                v2 = currvadj[currvidx[j]];
                v1 = currvadj[currvidx[j] + 1];
              }
              f2cnt = -1;
              v2chk = 1;
              for ( k = currvidx[v1 - 1]; k < currvidx[v1]; k++ )
              {
                if ( currvadj[k] == j + 1 )
                {
                  f2cnt = f2cnt + 1;
                  f2[f2cnt] = k;
                }
                if ( currvadj[k] == v2 )
                {
                  v2chk = 0;
                }
              }
              if ( v2chk == 1 )
              {
                for ( k = 0; k < f2cnt + 1; k++ )
                {
                  currvadj[f2[k]] = v2;
                  e1 = curreadj[f2[k]];
                }
              }
              else
              {
                t = 0;
                continue;
              }
              f2cnt = 0;
              for ( k = currvidx[v2 - 1]; k < currvidx[v2]; k++ )
              {
                if ( currvadj[k] == j + 1 )
                {
                  currvadj[k] = v1;
                  e2 = curreadj[k];
                  curreadj[k] = e1;
                  break;
                }
              }

              for ( k = currvidx[j]; k < currvidx[relcnt] - 2; k++ )
              {
                currvadj[k] = currvadj[k + 2];
                curreadj[k] = curreadj[k + 2];
              }
              for ( k = j + 2; k <= relcnt; k++ )
              {
                currvidx[k - 1] = currvidx[k] - 2;
              }
              for ( k = j + 1; k < relcnt; k++ )
              {
                currdeg[k - 1] = currdeg[k];
              }
              for ( k = e2; k <= Mcnt - 1; k++ )
              {
                einvmaptwo[k - 1] = einvmaptwo[k];
              }
              for ( k = j; k <= relcnt - 2; k++ )
              {
                vinvmaptwo[k] = vinvmaptwo[k + 1];
              }
              for ( k = 1; k <= origN; k++ )
              {
                if ( vmaptwo[k - 1] == j + 1 )
                {
                  vmaptwo[k - 1] = v2;
                }
              }
              for ( k = 1; k <= origN; k++ )
              {
                if ( vmaptwo[k - 1] > j + 1 )
                {
                  vmaptwo[k - 1]--;
                }
              }
              for ( k = 0; k < currvidx[relcnt - 1]; k++ )
              {
                if ( currvadj[k] > j + 1 )
                {
                  currvadj[k] = currvadj[k] - 1;
                }
                if ( curreadj[k] > e2 )
                {
                  curreadj[k] = curreadj[k] - 1;
                }
              }
              relcnt = relcnt - 1;
              Mcnt = Mcnt - 1;
              break;
            }
          }
        }
        if ( relcnt < 5 )
        {
          free( vmaptwo );
          free( currvadj );
          free( curreadj );
          free( currvidx );
          free( currdeg );
          free( f2 );
          free( einvmaptwo );
          free( vinvmaptwo );
          continue;
        }

        for ( j = 0; j < relcnt; j++ )
        {
          if ( currdeg[j] > maxdeg )
          {
            maxdeg = currdeg[j];
          }
        }
        int* Sixidx;
        Sixidx = (int*)calloc( 6 * Mcnt, sizeof( int ) );
        for ( j = 0; j < relcnt; j++ )
        {
          for ( k = currvidx[j]; k < currvidx[j + 1]; k++ )
          {
            if ( j + 1 < currvadj[k] )
            {
              Sixidx[curreadj[k] - 1] = j + 1;
              Sixidx[Mcnt + curreadj[k] - 1] = currvadj[k];
            }
            else
            {
              Sixidx[curreadj[k] - 1] = currvadj[k];
              Sixidx[Mcnt + curreadj[k] - 1] = j + 1;
            }
          }
        }

        int* incidente;
        incidente = (int*)calloc( relcnt, sizeof( int ) );
        for ( j = 0; j < Mcnt; j++ )
        {
          if ( incidente[Sixidx[j] - 1] == 0 )
          {
            incidente[Sixidx[j] - 1] = j + 1;
          }
          if ( incidente[Sixidx[Mcnt + j] - 1] == 0 )
          {
            incidente[Sixidx[Mcnt + j] - 1] = j + 1;
          }
        }

        int orig_relcnt = relcnt;
        int M_before_embedding = Mcnt;
        if ( embed_scheme == 1 )
        {
          long double* currx = (long double*)calloc( relcnt, sizeof( long double ) );
          long double* curry = (long double*)calloc( relcnt, sizeof( long double ) );
          Get_Spring_Data( Sixidx, currdeg, currvadj, currvidx, relcnt, Mcnt, &Clabel, &Dlabel, Cidx, relcnt, maxdeg, seed, 0, currx, curry );
          current_crossings = Cidx[Mcnt] / 2;
          free( currx );
          free( curry );
        }
        if ( embed_scheme == 2 )
        {
          Embed_Circ( &Clabel, &Dlabel, Cidx, Sixidx, relcnt, Mcnt, maxdeg, currvadj, curreadj, currvidx, seed );
          current_crossings = Cidx[Mcnt] / 2;
        }
        if ( embed_scheme == 3 )
        {
          int outputs[3];
          Embed_Planar_Subgraph( &Sixidx, &incidente, &Clabel, &Dlabel, &Cidx, relcnt, Mcnt, &currdeg, currvadj, currvidx, curreadj, seed, maxdeg, outputs, verbose );
          if ( Sixidx[2 * Mcnt] == 0 )
          {
            free( vmaptwo );
            free( currvadj );
            free( curreadj );
            free( currvidx );
            free( currdeg );
            free( f2 );
            free( Clabel );
            free( Dlabel );
            free( Cidx );
            free( Sixidx );
            free( incidente );
            free( einvmaptwo );
            free( vinvmaptwo );
            continue;
          }
          relcnt = outputs[0];
          Mcnt = outputs[1];
          current_crossings = outputs[2];
        }
        if ( embed_scheme == 5 )
        {
          // then we want the coords of the inverse mapping of the vertices.
          long double* currx = (long double*)calloc( relcnt, sizeof( long double ) );
          long double* curry = (long double*)calloc( relcnt, sizeof( long double ) );
          for ( j = 1; j <= relcnt; j++ )
          {
            currx[j - 1] = x[vinvmaptwo[j - 1] - 1];
            curry[j - 1] = y[vinvmaptwo[j - 1] - 1];
          }
          Get_Spring_Data( Sixidx, currdeg, currvadj, currvidx, relcnt, Mcnt, &Clabel, &Dlabel, Cidx, relcnt, maxdeg, seed, 1, currx, curry );
          current_crossings = Cidx[Mcnt] / 2;
          free( currx );
          free( curry );
        }
        acyclic = false;
        Quick_Cross_Main_Loop( &incidente, &Sixidx, &Clabel, &Dlabel, &Cidx, relcnt, Mcnt, &currdeg, current_crossings, min_scheme, bigface_depth, stop_crossings, true, orig_relcnt, outputs2, verbose );
        if ( outputs2[1] != M_before_embedding )
        {
          fprintf( stdout, "undosubd - not returned to original M \n" );
        }
        Mcnt = outputs2[1];
        relcnt = outputs2[2];
        Final_crossings = Final_crossings + outputs2[0];
        if ( fcnt + Cidx[Mcnt] > origM )
        {
          finClab1 = (int*)realloc( finClab1, ( fcnt + Cidx[Mcnt] ) * sizeof( int ) );
          finCidx1 = (int*)realloc( finCidx1, ( fcnt + Cidx[Mcnt] ) * sizeof( int ) );
        }

        for ( j = 0; j < Mcnt; j++ )
        {
          if ( vinvmaptwo[Sixidx[j] - 1] == el[einvmaptwo[j] - 1] || vinvmaptwo[Sixidx[Mcnt + j] - 1] == el[origM + einvmaptwo[j] - 1] )
          {
            for ( k = Cidx[j]; k < Cidx[j + 1]; k++ )
            {
              fcnt = fcnt + 1;
              finClab1[fcnt - 1] = einvmaptwo[Clabel[k] - 1];
              finCidx1[fcnt - 1] = einvmaptwo[j];
              finCcnt[einvmaptwo[j] - 1] = finCcnt[einvmaptwo[j] - 1] + 1;
            }
          }
          else if ( vinvmaptwo[Sixidx[j] - 1] == el[origM + einvmaptwo[j] - 1] || vinvmaptwo[Sixidx[Mcnt + j] - 1] == el[einvmaptwo[j] - 1] )
          {
            for ( k = Cidx[j + 1] - 1; k >= Cidx[j]; k-- )
            {
              fcnt = fcnt + 1;
              finClab1[fcnt - 1] = einvmaptwo[Clabel[k] - 1];
              finCidx1[fcnt - 1] = einvmaptwo[j];
              finCcnt[einvmaptwo[j] - 1] = finCcnt[einvmaptwo[j] - 1] + 1;
            }
          }
          else
          {
            fprintf( stdout, "map error 1" );
            getchar();
          }
        }
        free( vmaptwo );
        free( currvadj );
        free( curreadj );
        free( currvidx );
        free( currdeg );
        free( f2 );
        free( Clabel );
        free( Dlabel );
        free( Cidx );
        free( Sixidx );
        free( incidente );
        free( einvmaptwo );
        free( vinvmaptwo );
      }
    }
    else
    {
      if ( ctcnt > 1 )
      { // then more than one component
        int newN = ccidx[z] - ccidx[z - 1];
        int* relchk;
        relchk = (int*)malloc( newN * sizeof( int ) );
        int relchkmax = 0;
        for ( j = 0; j < newN; j++ )
        {
          relchk[j] = clst[ccidx[z - 1] + j];
          if ( relchk[j] > relchkmax )
          {
            relchkmax = relchk[j];
          }
        }
        int* vmaptwo; // this is vmaptwo in matlab.
        vmaptwo = (int*)calloc( origN, sizeof( int ) );
        int* currdeg;
        currdeg = (int*)malloc( newN * sizeof( int ) );
        int newM = 0;
        for ( j = 0; j < newN; j++ )
        {
          currdeg[j] = deg[relchk[j] - 1];
          newM = newM + currdeg[j];
        }
        newM = newM / 2;
        int* currvadj;
        int* curreadj;
        int* currSixidx;
        int* currvidx;
        currvadj = (int*)malloc( 2 * newM * sizeof( int ) );
        curreadj = (int*)malloc( 2 * newM * sizeof( int ) );
        currSixidx = (int*)calloc( 6 * newM, sizeof( int ) );
        currvidx = (int*)malloc( ( newN + 1 ) * sizeof( int ) );
        currvidx[0] = 0;
        for ( j = 0; j < newN; j++ )
        {

          currvidx[j + 1] = currvidx[j] + deg[relchk[j] - 1];
          vmaptwo[relchk[j] - 1] = j + 1;
        }
        int* einvmaptwo;
        einvmaptwo = (int*)calloc( newM, sizeof( int ) );
        int* echeck;
        echeck = (int*)calloc( M, sizeof( int ) );
        int Mcnt = 0;
        int kkcnt;
        for ( j = 0; j < newN; j++ )
        {
          kkcnt = -1;
          for ( k = vidx[relchk[j] - 1]; k < vidx[relchk[j]]; k++ )
          {
            if ( echeck[eadj[k] - 1] == 0 )
            {
              Mcnt = Mcnt + 1;
              echeck[eadj[k] - 1] = Mcnt;
              einvmaptwo[Mcnt - 1] = eadj[k];
            }
          }
          for ( k = currvidx[j]; k < currvidx[j + 1]; k++ )
          {
            kkcnt = kkcnt + 1;
            currvadj[k] = vmaptwo[vadj[vidx[clst[ccidx[z - 1] + j] - 1] + kkcnt] - 1];
            curreadj[k] = echeck[eadj[vidx[relchk[j] - 1] + kkcnt] - 1];
          }
        }
        int* Clabel;
        int* Dlabel;
        int* Cidx;
        Cidx = (int*)calloc( newM + 1, sizeof( int ) );
        int maxdeg = 0;
        Clabel = (int*)calloc( newN, sizeof( int ) );
        Dlabel = (int*)calloc( newN, sizeof( int ) );

        // start degree 2 reductions
        int t = 1;
        int v2;
        int f2cnt;
        int v2chk;
        int* f2;
        int e1;
        int e2;
        f2 = (int*)malloc( newN * sizeof( int ) );
        while ( t == 1 )
        {
          t = 0;
          if ( newN < 5 )
          {
            break;
          }
          for ( j = 0; j < newN; j++ )
          {
            if ( currdeg[j] == 2 )
            {
              t = 1;
              v1 = currvadj[currvidx[j]];
              v2 = currvadj[currvidx[j] + 1];
              if ( currdeg[v2 - 1] > 2 )
              {
                v2 = currvadj[currvidx[j]];
                v1 = currvadj[currvidx[j] + 1];
              }
              f2cnt = -1;
              v2chk = 1;
              for ( k = currvidx[v1 - 1]; k < currvidx[v1]; k++ )
              {
                if ( currvadj[k] == j + 1 )
                {
                  f2cnt = f2cnt + 1;
                  f2[f2cnt] = k;
                }
                if ( currvadj[k] == v2 )
                {
                  v2chk = 0;
                }
              }
              if ( v2chk == 1 )
              {
                for ( k = 0; k < f2cnt + 1; k++ )
                {
                  currvadj[f2[k]] = v2;
                  e1 = curreadj[f2[k]];
                }
              }
              else
              {
                t = 0;
                continue;
              }

              f2cnt = 0;
              for ( k = currvidx[v2 - 1]; k < currvidx[v2]; k++ )
              {
                if ( currvadj[k] == j + 1 )
                {
                  currvadj[k] = v1;
                  e2 = curreadj[k];
                  curreadj[k] = e1;
                  break;
                }
              }

              for ( k = currvidx[j]; k < currvidx[newN] - 2; k++ )
              {
                currvadj[k] = currvadj[k + 2];
                curreadj[k] = curreadj[k + 2];
              }

              for ( k = j + 2; k <= newN; k++ )
              {
                currvidx[k - 1] = currvidx[k] - 2;
              }

              for ( k = j + 1; k < newN; k++ )
              {
                currdeg[k - 1] = currdeg[k];
              }
              for ( k = e2; k <= newM - 1; k++ )
              {
                einvmaptwo[k - 1] = einvmaptwo[k];
              }
              for ( k = j; k <= newN - 2; k++ )
              {
                relchk[k] = relchk[k + 1];
              }
              for ( k = 1; k <= origN; k++ )
              {
                if ( vmaptwo[k - 1] == j + 1 )
                {
                  vmaptwo[k - 1] = v2;
                }
              }
              for ( k = 1; k <= origN; k++ )
              {
                if ( vmaptwo[k - 1] > j + 1 )
                {
                  vmaptwo[k - 1]--;
                }
              }
              for ( k = 0; k < currvidx[newN - 1]; k++ )
              {
                if ( currvadj[k] > j + 1 )
                {
                  currvadj[k] = currvadj[k] - 1;
                }
                if ( curreadj[k] > e2 )
                {
                  curreadj[k] = curreadj[k] - 1;
                }
              }

              newN = newN - 1;
              newM = newM - 1;
              break;
            }
          }
        }
        if ( newN < 5 )
        {
          free( visited );
          free( relchk );
          free( vmaptwo );
          free( currdeg );
          free( currvadj );
          free( curreadj );
          free( currSixidx );
          free( currvidx );
          free( einvmaptwo );
          free( echeck );
          free( f2 );
          continue;
        }
        for ( j = 0; j < newN; j++ )
        {
          for ( k = currvidx[j]; k < currvidx[j + 1]; k++ )
          {
            if ( j + 1 < currvadj[k] )
            {
              currSixidx[curreadj[k] - 1] = j + 1;
              currSixidx[newM + curreadj[k] - 1] = currvadj[k];
            }
            else
            {
              currSixidx[curreadj[k] - 1] = currvadj[k];
              currSixidx[newM + curreadj[k] - 1] = j + 1;
            }
          }
        }

        int* incidente;
        incidente = (int*)calloc( newN, sizeof( int ) );
        for ( j = 0; j < newM; j++ )
        {
          if ( incidente[currSixidx[j] - 1] == 0 )
          {
            incidente[currSixidx[j] - 1] = j + 1;
          }
          if ( incidente[currSixidx[newM + j] - 1] == 0 )
          {
            incidente[currSixidx[newM + j] - 1] = j + 1;
          }
        }

        for ( j = 0; j < newN; j++ )
        {
          if ( currdeg[j] > maxdeg )
          {
            maxdeg = currdeg[j];
          }
        }

        int orig_newN = newN;
        int M_before_embedding = newM;

        if ( embed_scheme == 1 )
        {
          long double* currx = (long double*)calloc( newN, sizeof( long double ) );
          long double* curry = (long double*)calloc( newN, sizeof( long double ) );
          Get_Spring_Data( currSixidx, currdeg, currvadj, currvidx, newN, newM, &Clabel, &Dlabel, Cidx, newN, maxdeg, seed, 0, currx, curry );
          current_crossings = Cidx[newM] / 2;
          free( currx );
          free( curry );
        }
        if ( embed_scheme == 2 )
        {
          Embed_Circ( &Clabel, &Dlabel, Cidx, currSixidx, newN, newM, maxdeg, currvadj, curreadj, currvidx, seed );
          current_crossings = Cidx[newM] / 2;
        }
        if ( embed_scheme == 3 )
        {
          int outputs[3];
          Embed_Planar_Subgraph( &currSixidx, &incidente, &Clabel, &Dlabel, &Cidx, newN, newM, &currdeg, currvadj, currvidx, curreadj, seed, maxdeg, outputs, verbose );
          if ( currSixidx[2 * newM] == 0 )
          {
            free( visited );
            free( relchk );
            free( vmaptwo );
            free( currdeg );
            free( currvadj );
            free( curreadj );
            free( currSixidx );
            free( currvidx );
            free( einvmaptwo );
            free( echeck );
            free( f2 );
            free( incidente );
            free( Clabel );
            free( Dlabel );
            free( Cidx );
            continue;
          }
          newN = outputs[0];
          newM = outputs[1];
          current_crossings = outputs[2];
        }
        if ( embed_scheme == 5 )
        {
          // then we want the coords of the inverse mapping of the vertices.
          long double* currx = (long double*)calloc( newN, sizeof( long double ) );
          long double* curry = (long double*)calloc( newN, sizeof( long double ) );
          for ( j = 1; j <= newN; j++ )
          {
            currx[j - 1] = x[vinvmapdel1[relchk[j - 1] - 1] - 1];
            curry[j - 1] = y[vinvmapdel1[relchk[j - 1] - 1] - 1];
          }
          Get_Spring_Data( currSixidx, currdeg, currvadj, currvidx, newN, newM, &Clabel, &Dlabel, Cidx, newN, maxdeg, seed, 1, currx, curry );
          current_crossings = Cidx[newM] / 2;
          free( currx );
          free( curry );
        }
        acyclic = false;
        Quick_Cross_Main_Loop( &incidente, &currSixidx, &Clabel, &Dlabel, &Cidx, newN, newM, &currdeg, current_crossings, min_scheme, bigface_depth, stop_crossings - Final_crossings, z != ctcnt, orig_newN, outputs2, verbose );
        if ( outputs2[1] != M_before_embedding )
        {
          fprintf( stdout, "undosubd2 - not returned to original M" );
        }
        newM = outputs2[1];
        newN = outputs2[2];
        Final_crossings = Final_crossings + outputs2[0];
        if ( fcnt + Cidx[newM] > origM )
        {
          finClab1 = (int*)realloc( finClab1, ( fcnt + Cidx[newM] ) * sizeof( int ) );
          finCidx1 = (int*)realloc( finCidx1, ( fcnt + Cidx[newM] ) * sizeof( int ) );
        }
        for ( j = 0; j < newM; j++ )
        {
          if ( vinvmapdel1[relchk[currSixidx[j] - 1] - 1] == el[einvmapdel1[einvmaptwo[j] - 1] - 1] || vinvmapdel1[relchk[currSixidx[newM + j] - 1] - 1] == el[origM + einvmapdel1[einvmaptwo[j] - 1] - 1] )
          {
            for ( k = Cidx[j]; k < Cidx[j + 1]; k++ )
            {
              fcnt = fcnt + 1;
              finClab1[fcnt - 1] = einvmaptwo[Clabel[k] - 1];
              finCidx1[fcnt - 1] = einvmaptwo[j];
              finCcnt[einvmaptwo[j] - 1] = finCcnt[einvmaptwo[j] - 1] + 1;
            }
          }
          else if ( vinvmapdel1[relchk[currSixidx[newM + j] - 1] - 1] == el[einvmapdel1[einvmaptwo[j] - 1] - 1] || vinvmapdel1[relchk[currSixidx[j] - 1] - 1] == el[origM + einvmapdel1[einvmaptwo[j] - 1] - 1] )
          {
            for ( k = Cidx[j + 1] - 1; k >= Cidx[j]; k-- )
            {
              fcnt = fcnt + 1;
              finClab1[fcnt - 1] = einvmaptwo[Clabel[k] - 1];
              finCidx1[fcnt - 1] = einvmaptwo[j];
              finCcnt[einvmaptwo[j] - 1] = finCcnt[einvmaptwo[j] - 1] + 1;
            }
          }
          else
          {
            fprintf( stdout, "map error 2" );
            getchar();
          }
        }
        free( relchk );
        free( vmaptwo );
        free( currdeg );
        free( currvadj );
        free( curreadj );
        free( currSixidx );
        free( currvidx );
        free( f2 );
        free( incidente );
        free( Clabel );
        free( Dlabel );
        free( Cidx );
        free( echeck );
        free( einvmaptwo );
      }
      else
      { // single component, so just use initial info.
        // deg 2 reductions
        int t = 1;
        int v2;
        int f2cnt;
        int v2chk;
        int* f2;
        int e1;
        int e2;
        f2 = (int*)malloc( N * sizeof( int ) );
        int* Clabel;
        int* Dlabel;
        int* Cidx;
        int maxdeg = 0;
        Clabel = (int*)calloc( N, sizeof( int ) );
        Dlabel = (int*)calloc( N, sizeof( int ) );
        Cidx = (int*)calloc( ( M + 1 ), sizeof( int ) );
        int* vtwo;
        vtwo = (int*)calloc( origN, sizeof( int ) );
        for ( j = 1; j <= N; j++ )
        {
          vtwo[vinvmapdel1[j - 1] - 1] = j;
        }

        while ( t == 1 )
        {
          t = 0;
          if ( N < 5 )
          {
            break;
          }
          for ( j = 0; j < N; j++ )
          {
            if ( deg[j] == 2 )
            {
              t = 1;
              v1 = vadj[vidx[j]];
              v2 = vadj[vidx[j] + 1];
              if ( deg[v2 - 1] > 2 )
              {
                v2 = vadj[vidx[j]];
                v1 = vadj[vidx[j] + 1];
              }
              f2cnt = -1;
              v2chk = 1;
              for ( k = vidx[v1 - 1]; k < vidx[v1]; k++ )
              {
                if ( vadj[k] == j + 1 )
                {
                  f2cnt = f2cnt + 1;
                  f2[f2cnt] = k;
                }
                if ( vadj[k] == v2 )
                {
                  v2chk = 0;
                }
              }
              if ( v2chk == 1 )
              {
                for ( k = 0; k < f2cnt + 1; k++ )
                {
                  vadj[f2[k]] = v2;
                  e1 = eadj[f2[k]];
                }
              }
              else
              {
                t = 0;
                continue;
              }
              f2cnt = 0;
              for ( k = vidx[v2 - 1]; k < vidx[v2]; k++ )
              {
                if ( vadj[k] == j + 1 )
                {
                  vadj[k] = v1;
                  e2 = eadj[k];
                  eadj[k] = e1;
                  break;
                }
              }
              for ( k = vidx[j]; k < vidx[N] - 2; k++ )
              {
                vadj[k] = vadj[k + 2];
                eadj[k] = eadj[k + 2];
              }
              for ( k = j + 2; k <= N; k++ )
              {
                vidx[k - 1] = vidx[k] - 2;
              }
              for ( k = j + 1; k < N; k++ )
              {
                deg[k - 1] = deg[k];
              }
              for ( k = e2; k <= M - 1; k++ )
              {
                einvmapdel1[k - 1] = einvmapdel1[k];
              }
              for ( k = j; k <= N - 2; k++ )
              {
                vinvmapdel1[k] = vinvmapdel1[k + 1];
              }
              for ( k = 1; k <= origN; k++ )
              {
                if ( vtwo[k - 1] == j + 1 )
                {
                  vtwo[k - 1] = v2;
                }
              }
              for ( k = 1; k <= origN; k++ )
              {
                if ( vtwo[k - 1] > j + 1 )
                {
                  vtwo[k - 1]--;
                }
              }
              for ( k = 0; k < vidx[N - 1]; k++ )
              {
                if ( vadj[k] > j + 1 )
                {
                  vadj[k] = vadj[k] - 1;
                }
                if ( eadj[k] > e2 )
                {
                  eadj[k] = eadj[k] - 1;
                }
              }
              N = N - 1;
              M = M - 1;
              break;
            }
          }
        }
        if ( N < 5 )
        {
          free( visited );
          free( LP );
          free( Tcheck );
          free( T );
          free( p );
          free( clst2 );
          free( ccidx2 );
          free( elst );
          free( eidx );
          free( f2 );
          continue;
        }
        int* Sixidx;
        Sixidx = (int*)calloc( 6 * M, sizeof( int ) );
        for ( j = 0; j < N; j++ )
        {
          for ( k = vidx[j]; k < vidx[j + 1]; k++ )
          {
            if ( j + 1 < vadj[k] )
            {
              Sixidx[eadj[k] - 1] = j + 1;
              Sixidx[M + eadj[k] - 1] = vadj[k];
            }
            else
            {
              Sixidx[eadj[k] - 1] = vadj[k];
              Sixidx[M + eadj[k] - 1] = j + 1;
            }
          }
        }

        for ( j = 0; j < N; j++ )
        {
          if ( deg[j] > maxdeg )
          {
            maxdeg = deg[j];
          }
        }
        int* incidente;
        incidente = (int*)calloc( N, sizeof( int ) );
        for ( j = 0; j < M; j++ )
        {
          if ( incidente[Sixidx[j] - 1] == 0 )
          {
            incidente[Sixidx[j] - 1] = j + 1;
          }
          if ( incidente[Sixidx[M + j] - 1] == 0 )
          {
            incidente[Sixidx[M + j] - 1] = j + 1;
          }
        }

        int orig_N = N;
        int M_before_embedding = M;
        if ( embed_scheme == 1 )
        {
          long double* currx = (long double*)calloc( N, sizeof( long double ) );
          long double* curry = (long double*)calloc( N, sizeof( long double ) );
          Get_Spring_Data( Sixidx, deg, vadj, vidx, N, M, &Clabel, &Dlabel, Cidx, N, maxdeg, seed, 0, currx, curry );
          current_crossings = Cidx[M] / 2;
          free( currx );
          free( curry );
        }
        if ( embed_scheme == 2 )
        {
          Embed_Circ( &Clabel, &Dlabel, Cidx, Sixidx, N, M, maxdeg, vadj, eadj, vidx, seed );
          current_crossings = Cidx[M] / 2;
        }
        if ( embed_scheme == 3 )
        {
          for ( j = 0; j < M; j++ )
          {
            Cidx[j + 1] = 0;
          }
          int outputs[3];
          Embed_Planar_Subgraph( &Sixidx, &incidente, &Clabel, &Dlabel, &Cidx, N, M, &deg, vadj, vidx, eadj, seed, maxdeg, outputs, verbose );
          if ( Sixidx[2 * M] == 0 )
          {
            free( Sixidx );
            free( f2 );
            free( Clabel );
            free( Dlabel );
            free( Cidx );
            free( incidente );
            free( vtwo );
            continue;
          }
          N = outputs[0];
          M = outputs[1];
          current_crossings = outputs[2];
        }
        if ( embed_scheme == 5 )
        {
          // then we want the coords of the inverse mapping of the vertices.
          long double* currx = (long double*)calloc( N, sizeof( long double ) );
          long double* curry = (long double*)calloc( N, sizeof( long double ) );
          for ( j = 1; j <= N; j++ )
          {
            currx[j - 1] = x[vinvmapdel1[j - 1] - 1];
            curry[j - 1] = y[vinvmapdel1[j - 1] - 1];
          }
          Get_Spring_Data( Sixidx, deg, vadj, vidx, N, M, &Clabel, &Dlabel, Cidx, N, maxdeg, seed, 1, currx, curry );
          current_crossings = Cidx[M] / 2;
          free( currx );
          free( curry );
        }
        acyclic = false;
        Quick_Cross_Main_Loop( &incidente, &Sixidx, &Clabel, &Dlabel, &Cidx, N, M, &deg, current_crossings, min_scheme, bigface_depth, stop_crossings, false, orig_N, outputs2, verbose );

        if ( outputs2[1] != M_before_embedding )
        {
          fprintf( stdout, "undosubd3 - not returned to original M: %d, %d \n", M, outputs2[1] );
        }
        Final_crossings = outputs2[0];
        M = outputs2[1];
        N = outputs2[2];

        if ( fcnt + Cidx[M] > origM )
        {
          finClab1 = (int*)realloc( finClab1, ( fcnt + Cidx[M] ) * sizeof( int ) );
          finCidx1 = (int*)realloc( finCidx1, ( fcnt + Cidx[M] ) * sizeof( int ) );
        }

        for ( j = 0; j < M; j++ )
        {
          if ( vtwo[el[einvmapdel1[j] - 1] - 1] == Sixidx[j] && ( vtwo[el[origM + einvmapdel1[j] - 1] - 1] == Sixidx[M + j] || Sixidx[M + j] > N ) )
          {
            for ( k = Cidx[j]; k < Cidx[j + 1]; k++ )
            {
              fcnt = fcnt + 1;
              finClab1[fcnt - 1] = einvmapdel1[Clabel[k] - 1];
              finCidx1[fcnt - 1] = einvmapdel1[j];
              finCcnt[einvmapdel1[j] - 1] = finCcnt[einvmapdel1[j] - 1] + 1;
            }
          }
          else
          {
            if ( vtwo[el[origM + einvmapdel1[j] - 1] - 1] != Sixidx[j] || ( vtwo[el[einvmapdel1[j] - 1] - 1] != Sixidx[M + j] && Sixidx[M + j] <= N ) )
            {
              fprintf( stdout, "map error 3\n" );
              getchar();
            }
            for ( k = Cidx[j + 1] - 1; k >= Cidx[j]; k-- )
            {
              fcnt = fcnt + 1;
              finClab1[fcnt - 1] = einvmapdel1[Clabel[k] - 1];
              finCidx1[fcnt - 1] = einvmapdel1[j];
              finCcnt[einvmapdel1[j] - 1] = finCcnt[einvmapdel1[j] - 1] + 1;
            }
          }
        }

        free( Sixidx );
        free( f2 );
        free( Clabel );
        free( Dlabel );
        free( Cidx );
        free( incidente );
        free( vtwo );
      }
    }
    free( visited );
    free( LP );
    free( Tcheck );
    free( T );
    free( p );
    free( clst2 );
    free( ccidx2 );
    free( elst );
    free( eidx );
  }
  int* fdn;
  if ( fcnt > origM )
  {
    *finClab = (int*)realloc( *finClab, fcnt * sizeof( int ) );
  }
  fdn = (int*)calloc( origM, sizeof( int ) );
  for ( i = 0; i < origM; i++ )
  {
    finCidx[i + 1] = finCidx[i] + finCcnt[i];
  }
  for ( i = 0; i < fcnt; i++ )
  {
    fdn[finCidx1[i] - 1] = fdn[finCidx1[i] - 1] + 1;
    ( *finClab )[finCidx[finCidx1[i] - 1] - 1 + fdn[finCidx1[i] - 1]] = finClab1[i];
  }

  free( deg );
  free( vadj );
  free( degchk );
  free( ccidx );
  free( visited );
  free( stack );
  free( clst );
  free( finClab1 );
  free( finCcnt );
  free( finCidx1 );
  free( eadj );
  free( vinvmapdel1 );
  free( einvmapdel1 );
  free( fdn );
  finoutputs[0] = Final_crossings;

  if ( acyclic && verbose > 0 )
    fprintf( stdout, "Graph is acyclic, hence crossing number is 0.\n\n" );
}

void Create_Graph6( char* graph6, int len, int** elist, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Create_Graph6                                                          */
  /*                                                                         */
  /*  Creates the graph specified by the graph6 code given with flag -g      */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;

  int N;
  int edges = 0;
  int total_edges = 100;
  int* el1 = (int*)malloc( 100 * sizeof( int ) );
  int* el2 = (int*)malloc( 100 * sizeof( int ) );
  int* edge_next = (int*)malloc( 100 * sizeof( int ) );
  int edge_first = -1;

  for ( i = 0; i < 100; i++ )
    edge_next[i] = -1;

  int index = 0;
  int x = graph6[index];

  if ( x == 126 )
  {
    int x1 = graph6[1];
    if ( x1 == 126 )
    {
      x1 = graph6[2];
      int x2 = graph6[3] - 63;
      int x3 = graph6[4] - 63;
      int x4 = graph6[5] - 63;
      int x5 = graph6[6] - 63;
      int x6 = graph6[7] - 63;
      index = 8;

      int nums[6];
      nums[0] = x1;
      nums[1] = x2;
      nums[2] = x3;
      nums[3] = x4;
      nums[4] = x5;
      nums[5] = x6;
      N = dec_bin_dec( nums, 6 );
    }
    else
    {
      int x2 = graph6[2];
      int x3 = graph6[3];
      index = 4;

      int nums[3];
      nums[0] = x1 - 63;
      nums[1] = x2 - 63;
      nums[2] = x3 - 63;
      N = dec_bin_dec( nums, 3 );
    }
  }
  else
  {
    N = x - 63;
    index = 1;
  }

  int nums_left = N * ( N - 1 ) / 2;
  if ( len - index != ( nums_left + 5 ) / 6 )
  {
    fprintf( stderr, "Graph6 code is invalid! Code: %s\n", graph6 );
    exit( 1 );
  }
  int row = 2;
  int col = 1;

  while ( nums_left > 0 )
  {
    x = graph6[index] - 63;
    index++;

    bool b[6];
    dec_to_bin( x, b );

    for ( i = 0; i < 6; i++ )
    {
      if ( b[i] )
      {
        el1[edges] = col;
        el2[edges] = row;

        int v = edge_first;
        int prev_v = -1;

        while ( v != -1 )
        {
          if ( el1[v] < col || ( el1[v] == col && el2[v] < row ) )
          {
            prev_v = v;
            v = edge_next[v];
          }
          else
            break;
        }

        if ( prev_v == -1 )
        {
          edge_next[edges] = edge_first;
          edge_first = edges;
        }
        else
        {
          edge_next[prev_v] = edges;
          edge_next[edges] = v;
        }

        edges++;
        if ( edges == total_edges )
        {
          total_edges = total_edges + 100;
          el1 = (int*)realloc( el1, total_edges * sizeof( int ) );
          el2 = (int*)realloc( el2, total_edges * sizeof( int ) );
          edge_next = (int*)realloc( edge_next, total_edges * sizeof( int ) );
          for ( j = total_edges - 100; j < total_edges; j++ )
            edge_next[j] = -1;
        }
      }
      col++;
      if ( col == row )
      {
        col = 1;
        row++;
      }
    }

    nums_left = nums_left - 6;
  }

  *elist = (int*)malloc( 2 * edges * sizeof( int ) );

  int v = edge_first;
  int count = 0;
  while ( v != -1 )
  {
    ( *elist )[count] = el1[v];
    ( *elist )[count + edges] = el2[v];
    count++;
    v = edge_next[v];
  }

  free( el1 );
  free( el2 );
  free( edge_next );

  outputs[0] = N;
  outputs[1] = edges;
}

void Read_Graph( char* filename, int** elist, int* outputs )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Read_Graph                                                             */
  /*                                                                         */
  /*  Reads a graph from a file in edgelist format.                          */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;

  FILE* ifp;
  ifp = fopen( filename, "r" );
  if ( ifp == NULL )
  {
    fprintf( stderr, "Can't open input file %s\n", filename );
    exit( 1 );
  }
  int N = 0;
  int M = 0;
  int maxN = 100;
  int maxM = 1000;
  int* raw = (int*)malloc( 1000 * sizeof( int ) );
  int* next = (int*)malloc( 1000 * sizeof( int ) );
  int* first = (int*)malloc( 100 * sizeof( int ) );
  int* entries = (int*)malloc( 100 * sizeof( int ) );

  for ( i = 0; i < 100; i++ )
  {
    first[i] = -1;
    entries[i] = 0;
  }
  while ( true )
  {

    int v1 = 0;
    int v2 = 0;
    int check = fscanf( ifp, "%d %d", &v1, &v2 );
    if ( check == EOF )
      break;
    if ( check != 2 )
    {
      fprintf( stdout, "Format of %s is wrong. Each line should only contain two vertices.\n", filename );
      exit( 1 );
    }
    if ( v1 > v2 )
    {
      int temp = v1;
      v1 = v2;
      v2 = temp;
    }

    if ( v2 > N )
    {
      N = v2;
      if ( N > maxN )
      {
        maxN = maxN + N;
        first = (int*)realloc( first, maxN * sizeof( int ) );
        entries = (int*)realloc( entries, maxN * sizeof( int ) );
        for ( i = maxN - N; i < maxN; i++ )
        {
          first[i] = -1;
          entries[i] = 0;
        }
      }
    }
    int e = first[v1 - 1];
    bool insert_new_edge = true;
    for ( i = 0; i < entries[v1 - 1]; i++ )
    {
      if ( raw[e] == v2 )
      {
        insert_new_edge = false;
        break;
      }
      if ( i != entries[v1 - 1] - 1 )
        e = next[e];
    }
    if ( insert_new_edge )
    {
      if ( e == -1 )
        first[v1 - 1] = M;

      else
        next[e] = M;

      entries[v1 - 1] = entries[v1 - 1] + 1;
      raw[M] = v2;
      M++;

      if ( M >= maxM )
      {
        maxM = maxM + 1000;
        raw = (int*)realloc( raw, maxM * sizeof( int ) );
        next = (int*)realloc( next, maxM * sizeof( int ) );
      }
    }
  }
  fclose( ifp );

  *elist = (int*)malloc( 2 * M * sizeof( int ) );
  int count = 0;
  for ( i = 1; i <= N; i++ )
  {
    if ( entries[i - 1] == 0 )
      continue;
    int* temp = (int*)malloc( entries[i - 1] * sizeof( int ) );
    int e = first[i - 1];
    for ( j = 1; j <= entries[i - 1]; j++ )
    {
      temp[j - 1] = raw[e];
      e = next[e];
    }
    qsort( temp, entries[i - 1], sizeof( int ), ascend_cmpfunc );

    for ( j = 1; j <= entries[i - 1]; j++ )
    {
      ( *elist )[count + j - 1] = i;
      ( *elist )[count + j - 1 + M] = temp[j - 1];
    }

    count = count + entries[i - 1];
  }

  free( raw );
  free( next );
  free( first );
  free( entries );

  outputs[0] = N;
  outputs[1] = M;
}

void Create_Graph( int type, int p1, int p2, int p3, int** elist, int* outputs, int verbose )
{

  /***************************************************************************/
  /*                                                                         */
  /*  Create_Graph                                                           */
  /*                                                                         */
  /*  Creates a graph from one of the predefined types.                      */
  /*                                                                         */
  /*  Type 1 - Complete graph                                                */
  /*  Type 2 - Complete bipartite graph                                      */
  /*  Type 3 - Generalized Petersen graph                                    */
  /*                                                                         */
  /***************************************************************************/

  int i;
  int j;
  int N;
  int M;

  if ( type == 1 ) // complete
  {
    if ( verbose > 0 )
      fprintf( stdout, "Making complete graph K%d: optimal crossing number should be %d.\n", p1, ( p1 / 2 ) * ( ( p1 - 1 ) / 2 ) * ( ( p1 - 2 ) / 2 ) * ( ( p1 - 3 ) / 2 ) / 4 );

    N = p1;
    M = N * ( N - 1 ) / 2;
    ( *elist ) = (int*)malloc( 2 * M * sizeof( int ) );
    int count = 0;
    for ( i = 1; i <= N; i++ )
      for ( j = i + 1; j <= N; j++ )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = j;
        count++;
      }
  }
  if ( type == 2 ) // complete bipartite
  {
    if ( verbose > 0 )
      fprintf( stdout, "Making complete bipartite graph K%d,%d: optimal crossing number should be %d.\n", p1, p2, ( p1 / 2 ) * ( ( p1 - 1 ) / 2 ) * ( p2 / 2 ) * ( ( p2 - 1 ) / 2 ) );

    N = p1 + p2;
    M = p1 * p2;
    ( *elist ) = (int*)malloc( 2 * M * sizeof( int ) );
    int count = 0;
    for ( i = 1; i <= p1; i++ )
      for ( j = 1; j <= p2; j++ )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = p1 + j;
        count++;
      }
  }
  if ( type == 3 ) // Generalized Petersen Graph
  {
    if ( verbose > 0 )
      fprintf( stdout, "Making generalized Petersen graph GP(%d,%d).\n", p1, p2 );

    N = 2 * p1;
    M = 3 * p1;
    ( *elist ) = (int*)malloc( 2 * M * sizeof( int ) );
    int count = 0;
    for ( i = 1; i <= p1; i++ )
    {
      if ( i != p1 )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = i + 1;
        count++;
      }
      if ( i == 1 )
      {
        ( *elist )[count] = 1;
        ( *elist )[count + M] = p1;
        count++;
      }
      ( *elist )[count] = i;
      ( *elist )[count + M] = i + p1;
      count++;
    }

    for ( i = p1 + 1; i <= 2 * p1; i++ )
    {
      if ( i + p2 <= N )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = i + p2;
        count++;
      }
      if ( i - p2 <= p1 )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = i - p2 + p1;
        count++;
      }
    }
  }
  if ( type == 4 ) // Complete Tripartite Graph
  {
    if ( verbose > 0 )
    {
      int bound = 0;
      int n[3];
      n[0] = p1;
      n[1] = p2;
      n[2] = p3;
      int i;
      int j;
      int k;
      for ( i = 0; i < 3; i++ )
        for ( j = 0; j < 3; j++ )
          for ( k = j + 1; k < 3; k++ )
          {
            if ( i == j || i == k )
              continue;
            bound = bound + ( n[j] / 2 ) * ( ( n[j] - 1 ) / 2 ) * ( n[k] / 2 ) * ( ( n[k] - 1 ) / 2 ) + ( n[i] / 2 ) * ( ( n[i] - 1 ) / 2 ) * ( ( n[j] * n[k] ) / 2 );
          }
      fprintf( stdout, "Making complete tripartite graph K%d,%d,%d: optimal crossing number should be %d.\n", p1, p2, p3, bound );
    }

    N = p1 + p2 + p3;
    M = p1 * p2 + p1 * p3 + p2 * p3;
    ( *elist ) = (int*)malloc( 2 * M * sizeof( int ) );

    int count = 0;
    for ( i = 1; i <= p1; i++ )
      for ( j = 1; j <= p2; j++ )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = p1 + j;
        count++;
      }

    for ( i = 1; i <= p1; i++ )
      for ( j = 1; j <= p3; j++ )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = p1 + p2 + j;
        count++;
      }

    for ( i = 1; i <= p2; i++ )
      for ( j = 1; j <= p3; j++ )
      {
        ( *elist )[count] = p1 + i;
        ( *elist )[count + M] = p1 + p2 + j;
        count++;
      }
  }
  if ( type == 5 )
  {
    if ( verbose > 0 )
      fprintf( stdout, "Making Sheehan graph S_%d.\n", p1 );

    N = p1;
    M = N * N / 4 + 1;
    ( *elist ) = (int*)malloc( 2 * M * sizeof( int ) );

    int count = 0;
    for ( i = 1; i < N; i++ )
    {
      if ( i % 2 == 1 )
      {
        ( *elist )[count] = i;
        ( *elist )[count + M] = i + 1;
        count++;

        if ( i == 1 )
        {
          ( *elist )[count] = i;
          ( *elist )[count + M] = N;
          count++;
        }
      }
      else
      {
        for ( j = i + 1; j <= N; j++ )
        {
          ( *elist )[count] = i;
          ( *elist )[count + M] = j;
          count++;
        }
      }
    }
  }
  if ( type == 6 )
  {
    if ( verbose > 0 )
    {
      int optc = p1;
      if ( p1 == 5 )
        optc = 4;
      fprintf( stdout, "Making Flower Snark I_%d: optimal crossing number should be %d.\n", p1, optc );
    }

    N = 4 * p1;
    M = 6 * p1;
    ( *elist ) = (int*)malloc( 2 * M * sizeof( int ) );

    int count = 0;
    for ( i = 1; i <= p1; i++ )
    {
      ( *elist )[count] = 4 * ( i - 1 ) + 1;
      ( *elist )[count + M] = 4 * ( i - 1 ) + 3;
      count++;

      ( *elist )[count] = 4 * ( i - 1 ) + 2;
      ( *elist )[count + M] = 4 * ( i - 1 ) + 3;
      count++;

      ( *elist )[count] = 4 * ( i - 1 ) + 3;
      ( *elist )[count + M] = 4 * ( i - 1 ) + 4;
      count++;
    }

    for ( i = 1; i <= p1 - 1; i++ )
    {
      ( *elist )[count] = 4 * i;
      ( *elist )[count + M] = 4 * i + 4;
      count++;
    }

    ( *elist )[count] = 4;
    ( *elist )[count + M] = 4 * p1;
    count++;

    for ( i = 1; i <= p1 - 1; i++ )
    {
      ( *elist )[count] = 4 * i - 3;
      ( *elist )[count + M] = 4 * i + 1;
      count++;
    }

    ( *elist )[count] = 2;
    ( *elist )[count + M] = 4 * p1 - 3;
    count++;

    for ( i = 1; i <= p1 - 1; i++ )
    {
      ( *elist )[count] = 4 * i - 2;
      ( *elist )[count + M] = 4 * i + 2;
      count++;
    }

    ( *elist )[count] = 1;
    ( *elist )[count + M] = 4 * p1 - 2;
    count++;
  }

  outputs[0] = N;
  outputs[1] = M;
}

void Print_Help()
{

  /***************************************************************************/
  /*                                                                         */
  /*  Print_Help                                                             */
  /*                                                                         */
  /*  Prints out the command line argument options that may be used.         */
  /*                                                                         */
  /***************************************************************************/

  fprintf( stdout, "Quick Cross v1.0\n\n" );
  fprintf( stdout, "Correct format:    QC [-bcefmosvw]\n" );
  fprintf( stdout, "Exactly one of -c, -g or -f must be specified, all other flags are optional.\n\n" );
  fprintf( stdout, "Example: QC -c 1 12 -w output\n" );
  fprintf( stdout, "   Runs complete graph K_12 with default options and writes output to \"output\".\n\n\n" );
  fprintf( stdout, "Required flags:\n\n" );
  fprintf( stdout, "-f <filename>            Read edgelist from file.\n\n" );
  fprintf( stdout, "   OR\n\n" );
  fprintf( stdout, "-c <type> <a> [<b> <c>]  Use graph creator.\n" );
  fprintf( stdout, "                            Input 1: Type\n" );
  fprintf( stdout, "                              1 = Complete Graph K_a\n" );
  fprintf( stdout, "                              2 = Complete Bipartite Graph K_a,b\n" );
  fprintf( stdout, "                              3 = Generalized Petersen Graph GP(a,b)\n" );
  fprintf( stdout, "                              4 = Complete Tripartite Graph K_a,b,c\n" );
  fprintf( stdout, "                              5 = Sheehan Maximally Uniquely Hamil. Graph S_n\n" );
  fprintf( stdout, "                              6 = Flower Snark I_n (n must be odd and >= 5)\n" );
  fprintf( stdout, "\n   OR\n\n" );
  fprintf( stdout, "-g <code/filename>       Create graph from graph6 code. Alternatively, reads a\n" );
  fprintf( stdout, "                         set of graph6 codes from a file, with one graph6\n" );
  fprintf( stdout, "                         code per line. The graphs are solved in order.\n\n" );
  fprintf( stdout, "\n\nOptional flags:\n\n" );
  fprintf( stdout, "-b <depth>    [Def: 0]   Bigface depth. Ignored unless Min Scheme = 1.\n" );
  fprintf( stdout, "                              Warning: Setting this larger than 3 will be slow!\n\n" );
  fprintf( stdout, "-e <scheme>   [Def: 3]   Embedding Scheme.\n" );
  fprintf( stdout, "                              1 = Kanada spring model drawing.\n" );
  fprintf( stdout, "                              2 = Embed on a circle.\n" );
  fprintf( stdout, "                              3 = Use Planar embedding.\n" );
  fprintf( stdout, "                              4 = Provide a crossing list file using -q.\n" );
  fprintf( stdout, "                              5 = Provide a vertex coordinates file using -q.\n\n" );
  fprintf( stdout, "-m <scheme>   [Def: 2]   Minimisation Scheme.\n" );
  fprintf( stdout, "                              1 = Use Bigface scheme.\n" );
  fprintf( stdout, "                              2 = Use improvements when found.\n" );
  fprintf( stdout, "                              3 = Search for max improvement each iteration.\n\n" );
  fprintf( stdout, "-o <scheme>   [Def: 1]   Output Scheme.\n" );
  fprintf( stdout, "                              1 = Only print out final crossings.\n" );
  fprintf( stdout, "                              2 = Also print out crossings list.\n" );
  fprintf( stdout, "                              3 = Also print out planar embedding.\n\n" );
  fprintf( stdout, "-q <filename>            Secondary file input for use in embed schemes 4 and 5.\n\n" );
  fprintf( stdout, "-r <num>      [Def: 0]   Sets random seed. If seed=0 no random numbers are used.\n\n" );
  fprintf( stdout, "-s <num>      [Def: 0]   Stop if number of crossings is equal or less than num.\n\n" );
  fprintf( stdout, "-u <num>      [Def: -1]  Upper Bound.\n" );
  fprintf( stdout, "                              Continue running graph (incrementing seed each time)\n" );
  fprintf( stdout, "                              until an embedding is found with no larger number of\n" );
  fprintf( stdout, "                              crossings than the upper bound. If upper bound is set\n" );
  fprintf( stdout, "                              to be -1, this is ignored and the graph will only be\n" );
  fprintf( stdout, "                              run a single time.\n\n" );
  fprintf( stdout, "-v <level>    [Def: 0]   Verbosity.\n" );
  fprintf( stdout, "                              0 = Only display final output.\n" );
  fprintf( stdout, "                              1 = Display crossings at each iteration only.\n" );
  fprintf( stdout, "                              2 = Display all information.\n\n" );
  fprintf( stdout, "-w <filename>            Output to file.\n" );
  fprintf( stdout, "                              If not specified, all output is printed to screen.\n\n" );
  fprintf( stdout, "-z <seed_print>          Re-run the graph indefinitely and record minimum solution.\n" );
  fprintf( stdout, "                              Must provide a seed with '-r', does not work with embed\n" );
  fprintf( stdout, "                              schemes 4 or 5. seed_print is a non-negative integer \n" );
  fprintf( stdout, "                              which specifies the frequency of print outs (only \n" );
  fprintf( stdout, "                              displayed, not written to file.) \n\n" );
  exit( 1 );
}

int qc_main( int argc, char* argv[] )
{

  int i;
  int j;

  int seed = 0;

  int crossings_check = -1;

  int embed_scheme = 3;
  int min_scheme = 2;
  int bigface_depth = 0;
  int type = -1;
  int p1 = -1;
  int p2 = -1;
  int p3 = -1;
  int stop_crossings = 0;
  int verbose = 0;
  int output_scheme = 1;

  int upper_bound = -1;
  int print_seed = 0;
  bool dont_stop = false;
  bool use_create_graph = false;
  bool use_edgelist_file = false;
  bool use_graph6 = false;
  bool use_graph6file = false;

  char* input_filename;
  char* output_filename;
  char* graph6_filename;
  char* graph6;
  char* secondary_filename;
  bool write_output_file = false;

  bool run_graph = false;

  if ( argc == 1 )
  {
    Print_Help();
  }
  int arg = 1;
  while ( arg < argc )
  {
    char* token = argv[arg];
    if ( token[0] != '-' )
    {
      fprintf( stdout, "Incorrect format\n" );
      Print_Help();
    }
    arg++;
    if ( arg == argc )
      Print_Help();
    char c = token[1];
    bool valid_token = false;
    if ( c == 'f' || c == 'F' )
    {
      run_graph = true;
      use_edgelist_file = true;
      valid_token = true;
      if ( use_create_graph || use_graph6 )
      {
        fprintf( stdout, "Cannot use more than one of  -c, -g and -f!\n" );
        return 0;
      }
      input_filename = argv[arg];
      arg++;
    }
    if ( c == 'g' || c == 'G' )
    {
      run_graph = true;
      use_graph6 = true;
      valid_token = true;
      if ( use_edgelist_file || use_create_graph )
      {
        fprintf( stdout, "Cannot use more than one of  -c, -g and -f!\n" );
        return 0;
      }
      graph6 = argv[arg];
      arg++;
    }
    if ( c == 'q' || c == 'Q' )
    {
      valid_token = true;
      run_graph = true;
      secondary_filename = argv[arg];
      arg++;
    }
    if ( c == 'c' || c == 'C' )
    {
      run_graph = true;
      use_create_graph = true;
      if ( use_edgelist_file || use_graph6 )
      {
        fprintf( stdout, "Cannot use more than one of  -c, -g and -f!\n" );
        return 0;
      }
      valid_token = true;
      char* t = argv[arg];
      char n = t[0];
      if ( n != '1' && n != '2' && n != '3' && n != '4' && n != '5' && n != '6' )
      {
        fprintf( stderr, "Type must be either 1, 2, 3, 4, 5 or 6!\n" );
        return 1;
      }
      type = atoi( argv[arg] );
      arg++;
      if ( arg == argc )
      {
        Print_Help();
      }
      p1 = atoi( argv[arg] );
      if ( p1 <= 0 )
      {
        fprintf( stderr, "Create graph arguments must be positive integers!\n" );
        return 1;
      }
      if ( type == 1 && p1 <= 2 )
      {
        fprintf( stderr, "The first argument must be at least 3 for Complete graphs!\n" );
        return 1;
      }
      if ( type == 5 && p1 <= 4 )
      {
        fprintf( stderr, "The first argument must be at least 5 for Sheehan graphs!\n" );
        return 1;
      }
      if ( type == 6 && ( p1 <= 4 || p1 % 2 == 0 ) )
      {
        fprintf( stderr, "The first argument must be odd and at least 5 for Flower Snarks!\n" );
        return 1;
      }
      arg++;
      if ( type == 2 || type == 3 || type == 4 )
      {
        p2 = atoi( argv[arg] );
        if ( p2 <= 0 )
        {
          fprintf( stderr, "Create graph arguments must be positive integers!\n" );
          return 1;
        }
        if ( type == 3 && 2 * p2 >= p1 )
        {
          fprintf( stderr, "The second argument must be smaller than half the first for GP graphs!\n" );
          return 1;
        }
        if ( type == 3 && p1 <= 2 )
        {
          fprintf( stderr, "The first argument must be at least 3 for GP graphs!\n" );
          return 1;
        }
        arg++;
      }
      if ( type == 4 )
      {
        p3 = atoi( argv[arg] );
        if ( p3 <= 0 )
        {
          fprintf( stderr, "Create graph arguments must be positive integers!\n" );
          return 1;
        }
        arg++;
      }
    }
    if ( c == 'm' || c == 'M' )
    {
      valid_token = true;
      char* t = argv[arg];
      char n = t[0];
      if ( n != '1' && n != '2' && n != '3' )
      {
        fprintf( stderr, "Min Scheme must be either 1, 2 or 3!\n" );
        return 1;
      }

      min_scheme = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'e' || c == 'E' )
    {
      valid_token = true;
      char* t = argv[arg];
      char n = t[0];
      if ( n != '2' && n != '3' && n != '1' && n != '4' && n != '5' )
      {
        fprintf( stderr, "Embed Scheme must be 1,2,3 or 4\n" );
        return 1;
      }
      embed_scheme = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'b' || c == 'B' )
    {
      valid_token = true;
      bigface_depth = atoi( argv[arg] );
      if ( bigface_depth < 0 )
      {
        fprintf( stderr, "Bigface depth must be nonnegative!\n" );
        return 1;
      }
      arg++;
    }
    if ( c == 'v' || c == 'V' )
    {
      valid_token = true;
      char* t = argv[arg];
      char n = t[0];
      if ( n != '0' && n != '1' && n != '2' )
      {
        fprintf( stderr, "Verbose must be either 0, 1 or 2!\n" );
        return 1;
      }
      verbose = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'o' || c == 'O' )
    {
      valid_token = true;
      char* t = argv[arg];
      char n = t[0];
      if ( n != '1' && n != '2' && n != '3' )
      {
        fprintf( stderr, "Output Scheme must be 1, 2 or 3!\n" );
        return 1;
      }
      output_scheme = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'w' || c == 'W' )
    {
      valid_token = true;
      output_filename = argv[arg];
      write_output_file = true;
      arg++;
    }
    if ( c == 's' || c == 'S' )
    {
      valid_token = true;
      stop_crossings = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'r' || c == 'R' )
    {
      valid_token = true;
      seed = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'u' || c == 'U' )
    {
      valid_token = true;
      upper_bound = atoi( argv[arg] );
      arg++;
    }
    if ( c == 'z' || c == 'Z' )
    {
      valid_token = true;
      dont_stop = true;
      print_seed = atoi( argv[arg] );
      if ( print_seed < 1 )
      {
        fprintf( stderr, "print_seed must be nonnegative!\n" );
        return 1;
      }
      arg++;
    }
    if ( valid_token == false )
    {
      Print_Help();
    }
    // Options:
    // -f edgelist filename
    // -c create graph
    // -g graph6 format
    // -m min_scheme (2)
    // -e embed_scheme (3)
    // -b bigface_depth (0)
    // -v verbose (0)
    // -o output_scheme (1)
    // -w output file
    // -s stop_crossings (0)
    // -r random seed (0)
  }

  if ( upper_bound >= 0 && seed == 0 )
  {
    fprintf( stdout, "You must provide a nonnegative seed if you request an upper bound!" );
    return 1;
  }

  if ( upper_bound >= 0 && embed_scheme >= 4 )
  {
    fprintf( stdout, "Only embedding schemes 1, 2 or 3 may be used with the upper bound option!" );
    return 1;
  }

  if ( upper_bound >= 0 && stop_crossings > upper_bound )
  {
    fprintf( stdout, "Upper bound must be no larger than the stop crossings number!" );
    return 1;
  }
  if ( dont_stop == true )
  {
    if ( embed_scheme > 3 )
    {
      fprintf( stdout, "non-stop mode requires embed_scheme 1,2 or 3!" );
      return 1;
    }
    if ( seed == 0 )
    {
      fprintf( stdout, "non-stop mode requires a non-negative random seed!" );
      return 1;
    }
  }

  if ( use_graph6 )
  {
    FILE* g6f;
    if ( g6f = fopen( graph6, "r" ) )
    {
      use_graph6 = false;
      use_graph6file = true;
      graph6_filename = graph6;
      fclose( g6f );
    }
    else
    {
    }
  }
  int nonplancnt = 0;

  if ( use_graph6file )
  {
    FILE* g6f = fopen( graph6_filename, "r" );

    // FILE *crf = fopen(secondary_filename, "r");

    if ( write_output_file )
    {
      FILE* outf;
      outf = fopen( output_filename, "w" );
      fclose( outf );
    }

    char g6s[100000];
    bool firsttime = true;
    int count = 0;
    while ( fgets( g6s, 100000, g6f ) != NULL )
    {
      if ( seed > 0 )
        srand( seed );

      count++;
      int N;
      int M;
      int* elist;
      int Final_crossings2[1];
      Final_crossings2[0] = -1;

      int outputs[2];

      int g6s_len = strlen( g6s );
      char newg6s[g6s_len - 1];

      if ( g6s[strlen( g6s ) - 1] == '\n' )
      {
        for ( i = 0; i < strlen( g6s ) - 1; i++ )
        {
          newg6s[i] = g6s[i];
        }
        int lin;
        lin = 0;
        if ( g6s[strlen( g6s ) - 2] == '\r' ) // then we are in a linux system and end of line character is different.
          lin = 1;

        Create_Graph6( newg6s, strlen( g6s ) - 1 - lin, &elist, outputs );
      }
      else
      {
        for ( i = 0; i < strlen( g6s ); i++ )
        {
          newg6s[i] = g6s[i];
        }
        int lin;
        lin = 0;
        if ( g6s[strlen( g6s ) - 1] == '\r' ) // then we are in a linux system and end of line character is different.
          lin = 1;

        Create_Graph6( newg6s, strlen( g6s ) - lin, &elist, outputs );
      }
      N = outputs[0];
      M = outputs[1];
      // fprintf(stdout,"read graph \n");
      int* finClabel = (int*)calloc( M, sizeof( int ) );
      int* finCidx = (int*)calloc( M + 1, sizeof( int ) );
      long double* x;
      long double* y;
      int best_found = INT_MAX;
      int best_found_seed = 0;
      int rerun = 1;
      while ( rerun )
      {
        rerun = 0;

        if ( embed_scheme == 4 )
        {

          FILE* crf = fopen( secondary_filename, "r" );

          char crow[10000];
          if ( count == 1 )
          {
            crow[0] = '-';
          }
          int* origclab;
          int* origcidx;
          origclab = (int*)calloc( M, sizeof( int ) );
          origcidx = (int*)calloc( M + 1, sizeof( int ) );
          int currsize = N;
          int edgeid;
          int oldedgeid;
          char* val;
          int value;
          int freespot = 0;
          char delims1[] = ":";
          char delims2[] = " ";
          if ( crow[0] != 'C' )
          {
            fgets( crow, 10000, crf ); // first line is crossings line.
          }
          oldedgeid = 0;
          fgets( crow, 10000, crf ); // get first edge crossings line.

          while ( strcmp( crow, "\n" ) != 0 && crow[0] != 'C' )
          {
            val = strtok( crow, delims1 );
            edgeid = atoi( val ); // first value is the edge number, the rest are crossings on that edge.

            for ( i = oldedgeid + 2; i <= edgeid; i++ )
            {
              origcidx[i - 1] = freespot;
            }
            val = strtok( NULL, delims2 );
            value = atoi( val );
            while ( value != 0 )
            {

              freespot++;
              if ( freespot > currsize )
              {
                origclab = (int*)realloc( origclab, ( currsize + 100 ) * sizeof( int ) );
                currsize = currsize + 100;
              }
              origclab[freespot - 1] = value;
              val = strtok( NULL, delims2 );
              value = atoi( val );
            }
            origcidx[edgeid] = freespot;
            oldedgeid = edgeid;
            // get next line.
            if ( fgets( crow, 10000, crf ) == NULL )
            {
              break;
            }
          }

          for ( i = oldedgeid + 2; i <= M + 1; i++ )
          {
            origcidx[i - 1] = freespot;
          }

          fclose( crf );
          Final_crossings2[0] = Embedding_Given_Crossings( N, M, origclab, origcidx, elist, min_scheme, stop_crossings, bigface_depth, verbose, finClabel, finCidx );
          if ( Final_crossings2[0] == -1 )
          {
            return ( 1 );
          }
          free( origcidx );
          free( origclab );
        }
        else if ( embed_scheme == 5 )
        {
          FILE* crf = fopen( secondary_filename, "r" );
          x = (long double*)calloc( N, sizeof( long double ) );
          y = (long double*)calloc( N, sizeof( long double ) );
          char delims[] = " \t\r\n\v\f";
          char crow[10000];
          char* val;
          if ( count == 1 )
          {
            crow[0] = '-';
          }
          for ( i = 1; i <= N; i++ )
          {
            fgets( crow, 10000, crf );
            val = strtok( crow, delims );
            x[i - 1] = atof( val );
            val = strtok( NULL, delims );
            y[i - 1] = atof( val );
          }
          fgets( crow, 10000, crf ); // this should be an emtpy line.
          fclose( crf );

          Biconnected_Runner( elist, N, M, min_scheme, embed_scheme, bigface_depth, seed, verbose, output_scheme, stop_crossings, output_filename, Final_crossings2, &finClabel, finCidx, x, y );

          free( x );
          free( y );
        }
        else
        {
          Biconnected_Runner( elist, N, M, min_scheme, embed_scheme, bigface_depth, seed, verbose, output_scheme, stop_crossings, output_filename, Final_crossings2, &finClabel, finCidx, x, y );
        }

        if ( upper_bound != -1 && Final_crossings2[0] > upper_bound )
        {
          if ( verbose >= 1 )
            fprintf( stdout, "CROSSINGS LARGER THAN UPPER BOUND - INCREMENTING SEED TO %d AND STARTING AGAIN!\n\n", seed + 1 );
          rerun = 1;
          seed++;
          srand( seed );
        }
        if ( dont_stop == true )
        {
          if ( Final_crossings2[0] < best_found )
          {
            best_found = Final_crossings2[0];
            fprintf( stdout, "New best found: %d, at seed: %d\n", best_found, seed );
            best_found_seed = seed;
            FILE* outfile;

            if ( write_output_file )
            {
              outfile = fopen( output_filename, "w" );
              if ( firsttime )
                firsttime = false;
              else
                fprintf( outfile, "\n" );
              if ( outfile == NULL )
              {
                fprintf( stdout, "Error opening file! Printing output to stdout instead.\n" );
                write_output_file = false;
              }
            }
            else
            {
              if ( firsttime )
                firsttime = false;
              // else
              // fprintf(stdout,"\n");
            }

            if ( write_output_file )
              fprintf( outfile, "CR : %d\n", Final_crossings2[0] );
            else
              fprintf( stdout, "CR : %d\n", Final_crossings2[0] );

            if ( output_scheme == 2 )
            {
              for ( i = 0; i < M; i++ )
              {
                if ( finCidx[i + 1] - finCidx[i] != 0 )
                {
                  if ( write_output_file )
                    fprintf( outfile, "%d: ", i + 1 );
                  else
                    fprintf( stdout, "%d: ", i + 1 );

                  for ( j = finCidx[i]; j < finCidx[i + 1]; j++ )
                  {
                    if ( write_output_file )
                      fprintf( outfile, "%d", finClabel[j] );
                    else
                      fprintf( stdout, "%d", finClabel[j] );
                    if ( j != finCidx[i + 1] - 1 )
                      if ( write_output_file )
                        fprintf( outfile, " " );
                      else
                        fprintf( stdout, " " );
                  }
                  if ( write_output_file )
                    fprintf( outfile, "\n" );
                  else
                    fprintf( stdout, "\n" );
                }
              }
            }

            if ( output_scheme >= 3 )
            {
              int tt;
              int finN = N;
              int ce;
              int k;
              // then we compute a planarised graph and output that.
              int finM = M;
              int finM2 = M + finCidx[M];
              int* planvlabels = (int*)calloc( finCidx[M], sizeof( int ) );
              int* elistorig = (int*)calloc( finM2, sizeof( int ) );
              int* eidxorig = (int*)calloc( M + 1, sizeof( int ) );
              for ( i = 2; i <= M + 1; i++ )
              {
                eidxorig[i - 1] = i - 1 + finCidx[i - 1];
              }
              for ( i = 1; i <= M; i++ )
              {
                elistorig[eidxorig[i - 1]] = i;
              }
              int* finelist = (int*)calloc( 2 * finM2, sizeof( int ) );
              for ( i = 1; i <= M; i++ )
              {
                finelist[i - 1] = elist[i - 1];
                finelist[finM2 + i - 1] = elist[M + i - 1];
              }
              for ( i = 1; i <= M; i++ )
              {
                if ( finCidx[i] - finCidx[i - 1] > 0 )
                {
                  tt = finelist[finM2 + i - 1];
                }
                for ( j = 1; j <= finCidx[i] - finCidx[i - 1]; j++ )
                {
                  finM++;
                  elistorig[eidxorig[i - 1] + j] = finM;
                  if ( planvlabels[finCidx[i - 1] + j - 1] == 0 )
                  {
                    finN++;
                    planvlabels[finCidx[i - 1] + j - 1] = finN;
                    ce = finClabel[finCidx[i - 1] + j - 1];
                    for ( k = finCidx[ce - 1] + 1; k <= finCidx[ce]; k++ )
                    {
                      if ( finClabel[k - 1] == i )
                      {
                        planvlabels[k - 1] = finN;
                        break;
                      }
                    }
                  }
                }
                if ( finCidx[i] - finCidx[i - 1] > 0 )
                {
                  finelist[finM - 1] = tt;
                  finelist[finM + finM2 - 1] = planvlabels[finCidx[i] - 1];
                  finelist[i - 1 + finM2] = planvlabels[finCidx[i - 1]];
                }
                for ( j = 1; j <= finCidx[i] - finCidx[i - 1] - 1; j++ )
                {
                  if ( planvlabels[finCidx[i - 1] + j - 1] < planvlabels[finCidx[i - 1] + j] )
                  {
                    finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
                    finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
                  }
                  else
                  {
                    finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
                    finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
                  }
                }
              }
              if ( write_output_file )
              {
                for ( i = 1; i <= finM2; i++ )
                {
                  fprintf( outfile, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
                }
                fprintf( outfile, "\n" );
              }
              else
              {
                for ( i = 1; i <= finM2; i++ )
                {
                  fprintf( stdout, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
                }
                fprintf( stdout, "\n" );
              }
              free( planvlabels );
              free( elistorig );
              free( eidxorig );
              free( finelist );
            }
            if ( write_output_file )
            {
              fclose( outfile );
            }
          }
          if ( seed % print_seed == 0 )
          {
            fprintf( stdout, "Current best soln: %d, at seed: %d, current seed: %d \n", best_found, best_found_seed, seed );
          }
          rerun = 1;
          seed++;
          srand( seed );
        }
      }

      FILE* outfile;

      if ( write_output_file )
      {
        outfile = fopen( output_filename, "a" );
        if ( firsttime )
          firsttime = false;
        else
          fprintf( outfile, "\n" );
        if ( outfile == NULL )
        {
          fprintf( stdout, "Error opening file! Printing output to stdout instead.\n" );
          write_output_file = false;
        }
      }
      else
      {
        if ( firsttime )
          firsttime = false;
        // else
        // fprintf(stdout,"\n");
      }

      if ( write_output_file )
        fprintf( outfile, "CR : %d\n", Final_crossings2[0] );
      else
        fprintf( stdout, "CR : %d\n", Final_crossings2[0] );

      if ( output_scheme == 2 )
      {
        for ( i = 0; i < M; i++ )
        {
          if ( finCidx[i + 1] - finCidx[i] != 0 )
          {
            if ( write_output_file )
              fprintf( outfile, "%d: ", i + 1 );
            else
              fprintf( stdout, "%d: ", i + 1 );

            for ( j = finCidx[i]; j < finCidx[i + 1]; j++ )
            {
              if ( write_output_file )
                fprintf( outfile, "%d", finClabel[j] );
              else
                fprintf( stdout, "%d", finClabel[j] );
              if ( j != finCidx[i + 1] - 1 )
                if ( write_output_file )
                  fprintf( outfile, " " );
                else
                  fprintf( stdout, " " );
            }
            if ( write_output_file )
              fprintf( outfile, "\n" );
            else
              fprintf( stdout, "\n" );
          }
        }
      }

      if ( output_scheme >= 3 )
      {
        int tt;
        int finN = N;
        int ce;
        int k;
        // then we compute a planarised graph and output that.
        int finM = M;
        int finM2 = M + finCidx[M];
        int* planvlabels = (int*)calloc( finCidx[M], sizeof( int ) );
        int* elistorig = (int*)calloc( finM2, sizeof( int ) );
        int* eidxorig = (int*)calloc( M + 1, sizeof( int ) );
        for ( i = 2; i <= M + 1; i++ )
        {
          eidxorig[i - 1] = i - 1 + finCidx[i - 1];
        }
        for ( i = 1; i <= M; i++ )
        {
          elistorig[eidxorig[i - 1]] = i;
        }
        int* finelist = (int*)calloc( 2 * finM2, sizeof( int ) );
        for ( i = 1; i <= M; i++ )
        {
          finelist[i - 1] = elist[i - 1];
          finelist[finM2 + i - 1] = elist[M + i - 1];
        }
        for ( i = 1; i <= M; i++ )
        {
          if ( finCidx[i] - finCidx[i - 1] > 0 )
          {
            tt = finelist[finM2 + i - 1];
          }
          for ( j = 1; j <= finCidx[i] - finCidx[i - 1]; j++ )
          {
            finM++;
            elistorig[eidxorig[i - 1] + j] = finM;
            if ( planvlabels[finCidx[i - 1] + j - 1] == 0 )
            {
              finN++;
              planvlabels[finCidx[i - 1] + j - 1] = finN;
              ce = finClabel[finCidx[i - 1] + j - 1];
              for ( k = finCidx[ce - 1] + 1; k <= finCidx[ce]; k++ )
              {
                if ( finClabel[k - 1] == i )
                {
                  planvlabels[k - 1] = finN;
                  break;
                }
              }
            }
          }
          if ( finCidx[i] - finCidx[i - 1] > 0 )
          {
            finelist[finM - 1] = tt;
            finelist[finM + finM2 - 1] = planvlabels[finCidx[i] - 1];
            finelist[i - 1 + finM2] = planvlabels[finCidx[i - 1]];
          }
          for ( j = 1; j <= finCidx[i] - finCidx[i - 1] - 1; j++ )
          {
            if ( planvlabels[finCidx[i - 1] + j - 1] < planvlabels[finCidx[i - 1] + j] )
            {
              finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
              finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
            }
            else
            {
              finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
              finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
            }
          }
        }
        if ( write_output_file )
        {
          for ( i = 1; i <= finM2; i++ )
          {
            fprintf( outfile, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
          }
          fprintf( outfile, "\n" );
        }
        else
        {
          for ( i = 1; i <= finM2; i++ )
          {
            fprintf( stdout, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
          }
          fprintf( stdout, "\n" );
        }
        free( planvlabels );
        free( elistorig );
        free( eidxorig );
        free( finelist );
      }

      free( elist );
      free( finClabel );
      free( finCidx );

      if ( write_output_file )
        fclose( outfile );
    }

    fclose( g6f );
  }
  else
  {
    if ( run_graph )
    {
      if ( embed_scheme >= 4 )
      {
        fprintf( stdout, "Embed schemes 4 and 5 require graph6 format and secondary file input (-q option).\n" );
        return 1;
      }
      int N;
      int M;
      int* elist;
      int Final_crossings2[1];
      Final_crossings2[0] = -1;
      if ( use_create_graph )
      {
        int outputs[2];
        Create_Graph( type, p1, p2, p3, &elist, outputs, verbose );
        N = outputs[0];
        M = outputs[1];
      }
      if ( use_edgelist_file )
      {
        int outputs[2];
        Read_Graph( input_filename, &elist, outputs );
        N = outputs[0];
        M = outputs[1];
      }
      if ( use_graph6 )
      {
        int outputs[2];
        Create_Graph6( graph6, strlen( graph6 ), &elist, outputs );
        N = outputs[0];
        M = outputs[1];
      }
      int* finClabel = (int*)calloc( M, sizeof( int ) );
      int* finCidx = (int*)calloc( M + 1, sizeof( int ) );
      long double* x;
      long double* y;
      int* origclabel;
      int* origcidx;
      int best_found = INT_MAX;
      int best_found_seed = 0;
      if ( seed > 0 )
        srand( seed );

      int rerun = 1;
      while ( rerun )
      {
        rerun = 0;
        Biconnected_Runner( elist, N, M, min_scheme, embed_scheme, bigface_depth, seed, verbose, output_scheme, stop_crossings, output_filename, Final_crossings2, &finClabel, finCidx, x, y );
        if ( upper_bound >= 0 && Final_crossings2[0] > upper_bound )
        {
          if ( verbose >= 1 )
            fprintf( stdout, "CROSSINGS LARGER THAN UPPER BOUND - INCREMENTING SEED TO %d AND STARTING AGAIN!\n\n", seed + 1 );
          rerun = 1;
          seed++;
          srand( seed );
        }
        if ( dont_stop == true )
        {
          if ( Final_crossings2[0] < best_found )
          {
            best_found = Final_crossings2[0];
            fprintf( stdout, "New best found: %d, at seed: %d\n", best_found, seed );
            best_found_seed = seed;
            FILE* outfile;

            if ( write_output_file )
            {
              outfile = fopen( output_filename, "w" );
              if ( outfile == NULL )
              {
                fprintf( stdout, "Error opening file! Printing output to stdout instead.\n" );
                write_output_file = false;
              }
            }

            if ( write_output_file )
              fprintf( outfile, "CR : %d\n", Final_crossings2[0] );
            else
              fprintf( stdout, "CR : %d\n", Final_crossings2[0] );

            if ( output_scheme == 2 )
            {
              for ( i = 0; i < M; i++ )
              {
                if ( finCidx[i + 1] - finCidx[i] != 0 )
                {
                  if ( write_output_file )
                    fprintf( outfile, "%d: ", i + 1 );
                  else
                    fprintf( stdout, "%d: ", i + 1 );

                  for ( j = finCidx[i]; j < finCidx[i + 1]; j++ )
                  {
                    if ( write_output_file )
                      fprintf( outfile, "%d", finClabel[j] );
                    else
                      fprintf( stdout, "%d", finClabel[j] );
                    if ( j != finCidx[i + 1] - 1 )
                      if ( write_output_file )
                        fprintf( outfile, " " );
                      else
                        fprintf( stdout, " " );
                  }
                  if ( write_output_file )
                    fprintf( outfile, "\n" );
                  else
                    fprintf( stdout, "\n" );
                }
              }
            }

            if ( output_scheme >= 3 )
            {
              int tt;
              int finN = N;
              int ce;
              int k;
              // then we compute a planarised graph and output that.
              int finM = M;
              int finM2 = M + finCidx[M];
              int* planvlabels = (int*)calloc( finCidx[M], sizeof( int ) );
              int* elistorig = (int*)calloc( finM2, sizeof( int ) );
              int* eidxorig = (int*)calloc( M + 1, sizeof( int ) );
              for ( i = 2; i <= M + 1; i++ )
              {
                eidxorig[i - 1] = i - 1 + finCidx[i - 1];
              }
              for ( i = 1; i <= M; i++ )
              {
                elistorig[eidxorig[i - 1]] = i;
              }
              int* finelist = (int*)calloc( 2 * finM2, sizeof( int ) );
              for ( i = 1; i <= M; i++ )
              {
                finelist[i - 1] = elist[i - 1];
                finelist[finM2 + i - 1] = elist[M + i - 1];
              }
              for ( i = 1; i <= M; i++ )
              {
                if ( finCidx[i] - finCidx[i - 1] > 0 )
                {
                  tt = finelist[finM2 + i - 1];
                }
                for ( j = 1; j <= finCidx[i] - finCidx[i - 1]; j++ )
                {
                  finM++;
                  elistorig[eidxorig[i - 1] + j] = finM;
                  if ( planvlabels[finCidx[i - 1] + j - 1] == 0 )
                  {
                    finN++;
                    planvlabels[finCidx[i - 1] + j - 1] = finN;
                    ce = finClabel[finCidx[i - 1] + j - 1];
                    for ( k = finCidx[ce - 1] + 1; k <= finCidx[ce]; k++ )
                    {
                      if ( finClabel[k - 1] == i )
                      {
                        planvlabels[k - 1] = finN;
                        break;
                      }
                    }
                  }
                }
                if ( finCidx[i] - finCidx[i - 1] > 0 )
                {
                  finelist[finM - 1] = tt;
                  finelist[finM + finM2 - 1] = planvlabels[finCidx[i] - 1];
                  finelist[i - 1 + finM2] = planvlabels[finCidx[i - 1]];
                }
                for ( j = 1; j <= finCidx[i] - finCidx[i - 1] - 1; j++ )
                {
                  if ( planvlabels[finCidx[i - 1] + j - 1] < planvlabels[finCidx[i - 1] + j] )
                  {
                    finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
                    finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
                  }
                  else
                  {
                    finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
                    finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
                  }
                }
              }
              if ( write_output_file )
              {
                for ( i = 1; i <= finM2; i++ )
                {
                  fprintf( outfile, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
                }
                fprintf( outfile, "\n" );
              }
              else
              {
                for ( i = 1; i <= finM2; i++ )
                {
                  fprintf( stdout, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
                }
                fprintf( stdout, "\n" );
              }
              free( planvlabels );
              free( elistorig );
              free( eidxorig );
              free( finelist );
            }
            if ( write_output_file )
            {
              fclose( outfile );
            }
          }
          if ( seed % print_seed == 0 )
          {
            fprintf( stdout, "Current best soln: %d, at seed: %d, current seed: %d \n", best_found, best_found_seed, seed );
          }
          rerun = 1;
          seed++;
          srand( seed );
        }
      }

      FILE* outfile;

      if ( write_output_file )
      {
        outfile = fopen( output_filename, "w" );
        if ( outfile == NULL )
        {
          fprintf( stdout, "Error opening file! Printing output to stdout instead.\n" );
          write_output_file = false;
        }
      }

      if ( write_output_file )
        fprintf( outfile, "CR : %d\n", Final_crossings2[0] );
      else
        fprintf( stdout, "CR : %d\n", Final_crossings2[0] );

      if ( output_scheme == 2 )
      {
        for ( i = 0; i < M; i++ )
        {
          if ( finCidx[i + 1] - finCidx[i] != 0 )
          {
            if ( write_output_file )
              fprintf( outfile, "%d: ", i + 1 );
            else
              fprintf( stdout, "%d: ", i + 1 );

            for ( j = finCidx[i]; j < finCidx[i + 1]; j++ )
            {
              if ( write_output_file )
                fprintf( outfile, "%d", finClabel[j] );
              else
                fprintf( stdout, "%d", finClabel[j] );
              if ( j != finCidx[i + 1] - 1 )
                if ( write_output_file )
                  fprintf( outfile, " " );
                else
                  fprintf( stdout, " " );
            }
            if ( write_output_file )
              fprintf( outfile, "\n" );
            else
              fprintf( stdout, "\n" );
          }
        }
      }

      if ( output_scheme >= 3 )
      {
        int tt;
        int finN = N;
        int ce;
        int k;
        // then we compute a planarised graph and output that.
        int finM = M;
        int finM2 = M + finCidx[M];
        int* planvlabels = (int*)calloc( finCidx[M], sizeof( int ) );
        int* elistorig = (int*)calloc( finM2, sizeof( int ) );
        int* eidxorig = (int*)calloc( M + 1, sizeof( int ) );
        for ( i = 2; i <= M + 1; i++ )
        {
          eidxorig[i - 1] = i - 1 + finCidx[i - 1];
        }
        for ( i = 1; i <= M; i++ )
        {
          elistorig[eidxorig[i - 1]] = i;
        }
        int* finelist = (int*)calloc( 2 * finM2, sizeof( int ) );
        for ( i = 1; i <= M; i++ )
        {
          finelist[i - 1] = elist[i - 1];
          finelist[finM2 + i - 1] = elist[M + i - 1];
        }
        for ( i = 1; i <= M; i++ )
        {
          if ( finCidx[i] - finCidx[i - 1] > 0 )
          {
            tt = finelist[finM2 + i - 1];
          }
          for ( j = 1; j <= finCidx[i] - finCidx[i - 1]; j++ )
          {
            finM++;
            elistorig[eidxorig[i - 1] + j] = finM;
            if ( planvlabels[finCidx[i - 1] + j - 1] == 0 )
            {
              finN++;
              planvlabels[finCidx[i - 1] + j - 1] = finN;
              ce = finClabel[finCidx[i - 1] + j - 1];
              for ( k = finCidx[ce - 1] + 1; k <= finCidx[ce]; k++ )
              {
                if ( finClabel[k - 1] == i )
                {
                  planvlabels[k - 1] = finN;
                  break;
                }
              }
            }
          }
          if ( finCidx[i] - finCidx[i - 1] > 0 )
          {
            finelist[finM - 1] = tt;
            finelist[finM + finM2 - 1] = planvlabels[finCidx[i] - 1];
            finelist[i - 1 + finM2] = planvlabels[finCidx[i - 1]];
          }
          for ( j = 1; j <= finCidx[i] - finCidx[i - 1] - 1; j++ )
          {
            if ( planvlabels[finCidx[i - 1] + j - 1] < planvlabels[finCidx[i - 1] + j] )
            {
              finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
              finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
            }
            else
            {
              finelist[elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j];
              finelist[finM2 + elistorig[eidxorig[i - 1] + j] - 1] = planvlabels[finCidx[i - 1] + j - 1];
            }
          }
        }
        if ( write_output_file )
        {
          for ( i = 1; i <= finM2; i++ )
          {
            fprintf( outfile, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
          }
          fprintf( outfile, "\n" );
        }
        else
        {
          for ( i = 1; i <= finM2; i++ )
          {
            fprintf( stdout, "%d %d\n", finelist[i - 1], finelist[finM2 + i - 1] );
          }
          fprintf( stdout, "\n" );
        }
        free( planvlabels );
        free( elistorig );
        free( eidxorig );
        free( finelist );
      }

      if ( write_output_file )
        fclose( outfile );

      free( elist );
      free( finClabel );
      free( finCidx );
    }
  }
  return 0;
}
