#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define NRA 100 /* number of rows in matrix A */
#define NCA 100 /* number of columns in matrix A */
#define NCB 100 /* number of columns in matrix B */
#define MASTER 0 /* taskid of first task */
#define FROM_MASTER 1 /* setting a message type */
#define FROM_WORKER 10 /* setting a message type */

double** create_dynamic_arr(int rows, int cols) {
    double* mem = (double*)malloc(rows * cols * sizeof(double));
    double** A = (double**)malloc(rows * sizeof(double*));
    A[0] = mem;
    for (int i = 1; i < rows; i++) A[i] = A[i - 1] + cols;
    return A;
}

void free_dynamic_arr(double** arr) {
    free(arr[0]);
    free(arr);
}

int main(int argc, char* argv[]) {
    int numtasks,
            taskid,
            numworkers,
            i, j, rc = 0;
    double time;

    double** a = create_dynamic_arr(NRA, NCA);
    double** b = create_dynamic_arr(NCA, NCB);
    double** c = create_dynamic_arr(NRA, NCB);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    double** a_rows = create_dynamic_arr(NRA / numtasks, NRA);
    double** c_rows = create_dynamic_arr(NRA / numtasks, NRA);

    if (numtasks < 2) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }
    numworkers = numtasks - 1;

    if (taskid == MASTER) {
        printf("mpi_mm has started with %d tasks.\n", numtasks);

        for (i = 0; i < NRA; i++)
            for (j = 0; j < NCA; j++)
                a[i][j] = 10;
        for (i = 0; i < NCA; i++)
            for (j = 0; j < NCB; j++)
                b[i][j] = 10;

        time = MPI_Wtime();
    }

    MPI_Scatter(
            *a,
            NRA * NRA / numtasks,
            MPI_DOUBLE,
            *a_rows,
            NRA * NRA / numtasks,
            MPI_DOUBLE,
            MASTER,
            MPI_COMM_WORLD);

    MPI_Bcast(
            *b,
            NRA * NRA,
            MPI_DOUBLE,
            MASTER,
            MPI_COMM_WORLD);

    for (int i = 0; i < NRA / numtasks; i++) {
        for (int j = 0; j < NRA; j++) {
            c_rows[i][j] = 0.0e0;
        }
    }

    for (int i = 0; i < NRA / numtasks; i++) {
        for (int k = 0; k < NRA; k++) {
            for (int j = 0; j < NRA; j++) {
                c_rows[i][j] += a_rows[i][k] * b[k][j];
            }
        }
    }

    MPI_Gather(
            *c_rows,
            (NRA / numtasks) * NRA,
            MPI_DOUBLE,
            *c,
            (NRA / numtasks) * NRA,
            MPI_DOUBLE,
            MASTER,
            MPI_COMM_WORLD);

    if (taskid == MASTER) {
        time = MPI_Wtime() - time;

        printf("\nTotal time: %.2f\n", time);
    }

    free_dynamic_arr(a);
    free_dynamic_arr(b);
    free_dynamic_arr(c);
    free_dynamic_arr(a_rows);
    free_dynamic_arr(c_rows);

    MPI_Finalize();
}

//int main(int argc, char* argv[]) {
//    int numtasks, taskid, numworkers, source, dest, mtype, /* message type */
//            rows, /* rows of matrix A sent to each worker */
//            averow, extra, offset, /* used to determine rows sent to each worker*/
//            i, j, k, rc = 0;
//
//    double** a = create_dynamic_arr(NRA, NCA);
//    double** b = create_dynamic_arr(NCA, NCB);
//    double** c = create_dynamic_arr(NRA, NCB);
//
//    MPI_Status status;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
//    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
//
//    if (numtasks < 2) {
//        printf("Need at least two MPI tasks. Quitting...\n");
//        MPI_Abort(MPI_COMM_WORLD, rc);
//        exit(1);
//    }
//    numworkers = numtasks - 1;
//
//    if (taskid == MASTER) {
//        printf("mpi_mm has started with %d tasks.\n", numtasks);
//        for (i = 0; i < NRA; i++)
//            for (j = 0; j < NCA; j++)
//                a[i][j] = 10;
//        for (i = 0; i < NCA; i++)
//            for (j = 0; j < NCB; j++)
//                b[i][j] = 10;
//        double t1 = MPI_Wtime();
//        averow = NRA / numworkers;
//        extra = NRA % numworkers;
//        offset = 0;
//        for (dest = 1; dest <= numworkers; dest++) {
//            rows = (dest <= extra) ? averow + 1 : averow;
//            printf("Sending %d rows to task %d offset=%d\n", rows, dest, offset);
//            MPI_Send(&offset, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
//            MPI_Send(&rows, 1, MPI_INT, dest, FROM_MASTER + 1, MPI_COMM_WORLD);
//            MPI_Send(a[offset], rows * NCA, MPI_DOUBLE, dest, FROM_MASTER + 2,
//                     MPI_COMM_WORLD);
//            MPI_Send(b[0], NCA * NCB, MPI_DOUBLE, dest, FROM_MASTER + 3, MPI_COMM_WORLD);
//            offset = offset + rows;
//        }
//        /* Receive results from worker tasks */
//        for (source = 1; source <= numworkers; source++) {
//            MPI_Recv(&offset, 1, MPI_INT, source, FROM_WORKER, MPI_COMM_WORLD, &status);
//            MPI_Recv(&rows, 1, MPI_INT, source, FROM_WORKER + 1, MPI_COMM_WORLD, &status);
//            MPI_Recv(c[offset], rows * NCB, MPI_DOUBLE, source, FROM_WORKER + 2, MPI_COMM_WORLD, &status);
//            //printf("Received results from task %d\n", id);
//        }
//        t1 = MPI_Wtime() - t1;
//        printf("\nTotal time: %.2f\n", t1);
//    }
//        /******** worker task *****************/
//
//    else { /* if (taskid > MASTER) */
//        MPI_Recv(&offset, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
//        MPI_Recv(&rows, 1, MPI_INT, MASTER, FROM_MASTER + 1, MPI_COMM_WORLD, &status);
//        MPI_Recv(a[0], rows * NCA, MPI_DOUBLE, MASTER, FROM_MASTER + 2, MPI_COMM_WORLD, &status);
//        MPI_Recv(b[0], NCA * NCB, MPI_DOUBLE, MASTER, FROM_MASTER + 3, MPI_COMM_WORLD, &status);
//        for (k = 0; k < NCB; k++)
//            for (i = 0; i < rows; i++) {
//                c[i][k] = 0.0;
//                for (j = 0; j < NCA; j++)
//                    c[i][k] = c[i][k] + a[i][j] * b[j][k];
//            }
//        MPI_Send(&offset, 1, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
//        MPI_Send(&rows, 1, MPI_INT, MASTER, FROM_WORKER + 1, MPI_COMM_WORLD);
//        MPI_Send(c[0], rows * NCB, MPI_DOUBLE, MASTER, FROM_WORKER + 2, MPI_COMM_WORLD);
//    }
//
//    free_dynamic_arr(a);
//    free_dynamic_arr(b);
//    free_dynamic_arr(c);
//
//    MPI_Finalize();
//}

//int main(int argc, char* argv[]) {
//	int numtasks,
//		taskid,
//		numworkers,
//		source,
//		dest,
//		mtype, /* message type */
//		rows, /* rows of matrix A sent to each worker */
//		averow, extra, offset, /* used to determine rows sent to each worker*/
//		i, j, k, rc = 0;
//
//	double** a = create_dynamic_arr(NRA, NCA);
//	double** b = create_dynamic_arr(NCA, NCB);
//	double** c = create_dynamic_arr(NRA, NCB);
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
//	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
//	if (numtasks < 2) {
//		printf("Need at least two MPI tasks. Quitting...\n");
//		MPI_Abort(MPI_COMM_WORLD, rc);
//		exit(1);
//	}
//	numworkers = numtasks - 1;
//	if (taskid == MASTER) {
//		printf("mpi_mm has started with %d tasks.\n", numtasks);
//		for (i = 0; i < NRA; i++)
//			for (j = 0; j < NCA; j++)
//				a[i][j] = 10;
//		for (i = 0; i < NCA; i++)
//			for (j = 0; j < NCB; j++)
//				b[i][j] = 10;
//		double t1 = MPI_Wtime();
//		averow = NRA / numworkers;
//		extra = NRA % numworkers;
//		offset = 0;
//		MPI_Status recv_status1[NUM_PROCESSES * 2], recv_status2[NUM_PROCESSES];
//		MPI_Request send_req[NUM_PROCESSES * 4], recv_req1[NUM_PROCESSES * 2],
//			recv_req2[NUM_PROCESSES];
//		for (dest = 1; dest <= numworkers; dest++) {
//			rows = (dest <= extra) ? averow + 1 : averow;
//			printf("Sending %d rows to task %d offset=%d\n", rows, dest, offset);
//			MPI_Isend(&offset, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD, &send_req[(dest - 1) * 4]);
//			MPI_Isend(&rows, 1, MPI_INT, dest, FROM_MASTER + 1, MPI_COMM_WORLD, &send_req[(dest - 1) * 4 + 1]);
//			MPI_Isend(a[offset], rows * NCA, MPI_DOUBLE, dest, FROM_MASTER + 2,
//				MPI_COMM_WORLD, &send_req[(dest - 1) * 4 + 2]);
//			MPI_Isend(b[0], NCA * NCB, MPI_DOUBLE, dest, FROM_MASTER + 3, MPI_COMM_WORLD, &send_req[(dest - 1) * 4 + 3]);
//			offset = offset + rows;
//		}
//
//		int offsets_arr[NUM_PROCESSES];
//		int rows_arr[NUM_PROCESSES];
//
//		/* Receive results from worker tasks */
//		for (source = 1; source <= numworkers; source++) {
//			MPI_Irecv(&offsets_arr[source - 1], 1, MPI_INT, source, FROM_WORKER,
//				MPI_COMM_WORLD, &recv_req1[(source - 1) * 2]);
//			MPI_Irecv(&rows_arr[source - 1], 1, MPI_INT, source, FROM_WORKER + 1,
//				MPI_COMM_WORLD, &recv_req1[(source - 1) * 2 + 1]);
//		}
//
//		MPI_Waitall(numworkers * 2, &recv_req1[0], &recv_status1[0]);
//
//		for (source = 1; source <= numworkers; source++) {
//			MPI_Irecv(c[offsets_arr[source - 1]], rows_arr[source - 1] * NCB,
//				MPI_DOUBLE, source, FROM_WORKER + 2, MPI_COMM_WORLD,
//				&recv_req2[source - 1]);
//			printf("Received results from task %d\n", source);
//		}
//
//		MPI_Waitall(numworkers, &recv_req2[0], &recv_status2[0]);
//
//		t1 = MPI_Wtime() - t1;
//
//		printf("\nTotal time: %.2f\n", t1);
//	}
//	/******** worker task *****************/
//	else { /* if (taskid > MASTER) */
//		MPI_Status recv_status[4];
//		MPI_Request send_req[3], recv_req[4];
//		MPI_Irecv(&offset, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD,
//			&recv_req[0]);
//		MPI_Irecv(&rows, 1, MPI_INT, MASTER, FROM_MASTER + 1, MPI_COMM_WORLD,
//			&recv_req[1]);
//		MPI_Irecv(b[0], NCA * NCB, MPI_DOUBLE, MASTER, FROM_MASTER + 3,
//			MPI_COMM_WORLD, &recv_req[3]);
//		MPI_Wait(&recv_req[0], &recv_status[0]);
//		MPI_Wait(&recv_req[1], &recv_status[1]);
//
//		MPI_Isend(&offset, 1, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD,
//			&send_req[0]);
//		MPI_Isend(&rows, 1, MPI_INT, MASTER, FROM_WORKER + 1, MPI_COMM_WORLD,
//			&send_req[1]);
//
//		MPI_Irecv(a[0], rows * NCA, MPI_DOUBLE, MASTER, FROM_MASTER + 2,
//			MPI_COMM_WORLD, &recv_req[2]);
//
//		MPI_Wait(&recv_req[2], &recv_status[2]);
//		MPI_Wait(&recv_req[3], &recv_status[3]);
//
//		for (k = 0; k < NCB; k++)
//			for (i = 0; i < rows; i++) {
//				c[i][k] = 0.0;
//				for (j = 0; j < NCA; j++)
//					c[i][k] = c[i][k] + a[i][j] * b[j][k];
//			}
//
//		MPI_Isend(c[0], rows * NCB, MPI_DOUBLE, MASTER, FROM_WORKER + 2,
//			MPI_COMM_WORLD, &send_req[2]);
//
//		MPI_Status st;
//		MPI_Wait(&send_req[2], &st);
//	}
//
//	free_dynamic_arr(a);
//	free_dynamic_arr(b);
//	free_dynamic_arr(c);
//
//	MPI_Finalize();
//}