/*
   Running without arguments is equivalent to 1000 iterations with the
   5 celestial objects declared in the golden_bodies array.

   $ nbody.exe 1000 5

   The output of this shows the energy before and after the simulation,
   and should be:

   double:
   -0.169075164
   -0.169087605

   float:
   -0.169075206
   -0.169086471
   */

#include <cuda_runtime.h>
#include <cuda_occupancy.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <curand.h>
#include <curand_kernel.h>
#include <algorithm>

#define GPUTEST 1


using type = float;

const type pi{ 3.141592653589793 };
const type solar_mass{ 4 * pi * pi };
const type days_per_year{ 365.24 };

int blockSize;
int mingridSize;
int gridSize;

type *outData;

std::ofstream file;


template <typename T>
struct planet {
	T x, y, z;
	T vx, vy, vz;
	T mass;
};

struct planet<type> golden_bodies[5] = {
	{                               /* sun */
		0, 0, 0, 0, 0, 0, solar_mass
	},
	{                               /* jupiter */
		4.84143144246472090e+00,
		-1.16032004402742839e+00,
		-1.03622044471123109e-01,
		1.66007664274403694e-03 * days_per_year,
		7.69901118419740425e-03 * days_per_year,
		-6.90460016972063023e-05 * days_per_year,
		9.54791938424326609e-04 * solar_mass
	},
	{                               /* saturn */
		8.34336671824457987e+00,
		4.12479856412430479e+00,
		-4.03523417114321381e-01,
		-2.76742510726862411e-03 * days_per_year,
		4.99852801234917238e-03 * days_per_year,
		2.30417297573763929e-05 * days_per_year,
		2.85885980666130812e-04 * solar_mass
	},
	{                               /* uranus */
		1.28943695621391310e+01,
		-1.51111514016986312e+01,
		-2.23307578892655734e-01,
		2.96460137564761618e-03 * days_per_year,
		2.37847173959480950e-03 * days_per_year,
		-2.96589568540237556e-05 * days_per_year,
		4.36624404335156298e-05 * solar_mass
	},
	{                               /* neptune */
		1.53796971148509165e+01,
		-2.59193146099879641e+01,
		1.79258772950371181e-01,
		2.68067772490389322e-03 * days_per_year,
		1.62824170038242295e-03 * days_per_year,
		-9.51592254519715870e-05 * days_per_year,
		5.15138902046611451e-05 * solar_mass
	}
};

const type DT{ 1e-2 };
const type RECIP_DT{ static_cast<type>(1.0) / DT };


//velocity update for the kernals
template <typename T>
__global__ void adv_Update_GPU(int nbodies, planet<T> *bodies)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if (i < nbodies)
	{

		planet<T> &b1 = bodies[i];
		for (int j = i + 1; j < nbodies; j++)
		{
			planet<T> &b2 = bodies[j];
			T dx = b1.x - b2.x;
			T dy = b1.y - b2.y;
			T dz = b1.z - b2.z;
			T inv_distance = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
			T mag = inv_distance * inv_distance * inv_distance;
			b1.vx -= dx * b2.mass * mag;
			b1.vy -= dy * b2.mass * mag;
			b1.vz -= dz * b2.mass * mag;
			b2.vx += dx * b1.mass  * mag;
			b2.vy += dy * b1.mass  * mag;
			b2.vz += dz * b1.mass  * mag;
		}

		b1.x += b1.vx;
		b1.y += b1.vy;
		b1.z += b1.vz;
	}

	
}

template <typename T>
__global__ void scale_bodies_GPU(int nbodies, planet<T> *bodies, T scale)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if (i < nbodies)
	{
		bodies[i].mass *= scale*scale;
		bodies[i].vx *= scale;
		bodies[i].vy *= scale;
		bodies[i].vz *= scale;
	}
}

//template <typename T>
//__global__ void energy_GPU(int nbodies, T *addReduc, T *subReduc, planet<T> *bodies)
//{
//	extern __shared__ T e[];
//
//	//T e = 0.0;
//	unsigned int threadID = threadIdx.x;
//
//	unsigned int i = threadIdx.x + blockIdx.x*blockDim.x;
//
//
//
//	if (i < nbodies)
//	{
//		planet<T> &b = bodies[i];
//		e[threadID] = 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz);
//	}
//
//
//	for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1)
//	{
//		if (threadID < stride)
//		{
//			e[threadID] += e[threadID + stride];
//		}
//
//	}
//	if (threadID == 0)
//	{
//		addReduc[blockIdx.x] = e[0];
//	}
//
//
//	e[threadID] = 0;
//
//	if (i < nbodies)
//	{
//		for (int iter = i + 1; iter < nbodies; iter++){
//			planet<T> &b = bodies[i];
//			planet<T> &b2 = bodies[iter];
//			T dx = b.x - b2.x;
//			T dy = b.y - b2.y;
//			T dz = b.z - b2.z;
//			T distance = sqrt(dx * dx + dy * dy + dz * dz);
//			T var = ((b.mass * b2.mass) / distance);
//			e[threadID] += var;
//		}
//	}
//
//	__syncthreads();
//
//	for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1)
//	{
//		if (threadID < stride)
//		{
//			e[threadID] += e[threadID + stride];
//		}
//		__syncthreads();
//	}
//
//	if (threadID == 0)
//	{
//		subReduc[blockIdx.x] = e[0];
//	}
//}

template <typename T>
__global__ void offset_momentum_GPU(planet<T> *bodies, T *outData, int nbodies, int d_gridSize, T solarMass)
{

	T px = 0.0, py = 0.0, pz = 0.0;
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int doCount = 0;
	reduction_offset_mom(nbodies, bodies, outData, 0);

	if (doCount < 3)
	{
		if (i > 0 && i < d_gridSize)
		{
			outData[0] += outData[i];
		}
		if (doCount == 0)
		px = outData[0];
		else if (doCount == 1)
			py = outData[0];
		else if (doCount == 2)
			pz = outData[0];
		doCount++;
	}

	bodies[0].vx = -px / solarMass;
	bodies[0].vy = -py / solarMass;
	bodies[0].vz = -pz / solarMass;
}

template <typename T>
__device__ void reduction_offset_mom(int nbodies, planet<T> *bodies, T *outData, int step)
{
	extern __shared__ T sharedData[];
	// each thread loads one element from global to shared mem
	unsigned int threadID = threadIdx.x;
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (step == 0)
	{
		if (i < nbodies)
			sharedData[threadID] = bodies[i].vx * bodies[i].mass;

	}
	else if (step == 1)
	{
		if (i < nbodies)
			sharedData[threadID] = bodies[i].vy * bodies[i].mass;
	}
	else if (step == 2)
	{
		if (i < nbodies)
			sharedData[threadID] = bodies[i].vz * bodies[i].mass;
	}


	// do reduction in shared mem
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (threadID < s) {
			sharedData[threadID] += sharedData[threadID + s];
		}

	}
	// write result for this block to global mem
	if (threadID == 0)
		outData[blockIdx.x] = sharedData[0];
}


template <typename T>
void advance(int nbodies, planet<T> *bodies)
{
	int i, j;

	for (i = 0; i < nbodies; ++i) {
		planet<T> &b = bodies[i];
		for (j = i + 1; j < nbodies; j++) {
			planet<T> &b2 = bodies[j];
			T dx = b.x - b2.x;
			T dy = b.y - b2.y;
			T dz = b.z - b2.z;
			T inv_distance = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
			T mag = inv_distance * inv_distance * inv_distance;
			b.vx -= dx * b2.mass * mag;
			b.vy -= dy * b2.mass * mag;
			b.vz -= dz * b2.mass * mag;
			b2.vx += dx * b.mass  * mag;
			b2.vy += dy * b.mass  * mag;
			b2.vz += dz * b.mass  * mag;
		}
	}
	for (i = 0; i < nbodies; ++i) {
		planet<T> &b = bodies[i];
		b.x += b.vx;
		b.y += b.vy;
		b.z += b.vz;
	}
}

template <typename T>
T energy(int nbodies, planet<T> *bodies)
{
	T e = 0.0;
	for (int i = 0; i < nbodies; ++i) {
		planet<T> &b = bodies[i];
		e += 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz);
		for (int j = i + 1; j < nbodies; j++) {
			planet<T> &b2 = bodies[j];
			T dx = b.x - b2.x;
			T dy = b.y - b2.y;
			T dz = b.z - b2.z;
			T distance = sqrt(dx * dx + dy * dy + dz * dz);
			e -= (b.mass * b2.mass) / distance;
		}
	}
	return e;
}

template <typename T>
void offset_momentum(int nbodies, planet<T> *bodies)
{
	T px = 0.0, py = 0.0, pz = 0.0;
	for (int i = 0; i < nbodies; ++i) {
		px += bodies[i].vx * bodies[i].mass;
		py += bodies[i].vy * bodies[i].mass;
		pz += bodies[i].vz * bodies[i].mass;
	}
	bodies[0].vx = -px / solar_mass;
	bodies[0].vy = -py / solar_mass;
	bodies[0].vz = -pz / solar_mass;
}

/*
 * Rescale certain properties of bodies. That allows doing
 * consequential advance()'s as if dt were equal to 1.0.
 *
 * When all advances done, rescale bodies back to obtain correct energy.
 */
template <typename T>
void scale_bodies(int nbodies, planet<T> *bodies, T scale)
{
	for (int i = 0; i < nbodies; ++i) {
		bodies[i].mass *= scale*scale;
		bodies[i].vx *= scale;
		bodies[i].vy *= scale;
		bodies[i].vz *= scale;
	}
}

template <typename T>
void init_random_bodies(int nbodies, planet<T> *bodies)
{
	for (int i = 0; i < nbodies; ++i) {
		bodies[i].x = (T)rand() / RAND_MAX;
		bodies[i].y = (T)rand() / RAND_MAX;
		bodies[i].z = (T)rand() / RAND_MAX;
		bodies[i].vx = (T)rand() / RAND_MAX;
		bodies[i].vy = (T)rand() / RAND_MAX;
		bodies[i].vz = (T)rand() / RAND_MAX;
		bodies[i].mass = (T)rand() / RAND_MAX;
	}
}

int isPowerOfTwo(unsigned int x)
{
	return ((x != 0) && !(x & (x - 1)));
}

template <typename T>
void callOffSet(int nbodies, planet<T> *Gbodies)
{
	T *d_Array;

	cudaMalloc(&d_Array, gridSize*sizeof(T));

	offset_momentum_GPU << <gridSize, blockSize, nbodies * sizeof(T) >> >(Gbodies, d_Array, nbodies, gridSize, solar_mass); //need to fix outData

	cudaFree(d_Array);
}

//template <typename T>
//T callEnergy(int nbodies, planet<T> *Gbodies)
//{
//	
//	T *h_addArray = new T[gridSize];
//	T *h_subArray = new T[gridSize];
//
//	T *d_addArray; cudaMalloc((void**)&d_addArray, gridSize * sizeof(T));
//	T *d_subArray; cudaMalloc((void**)&d_subArray, gridSize * sizeof(T));
//
//	energy_GPU << <gridSize, blockSize, nbodies * sizeof(T) >> >(nbodies, d_addArray, d_subArray, Gbodies);
//	cudaMemcpy(h_addArray, d_addArray, gridSize * sizeof(T), cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_subArray, d_subArray, gridSize * sizeof(T), cudaMemcpyDeviceToHost);
//
//	for (int i = 1; i < gridSize; i++){
//		h_addArray[0] += h_addArray[i];
//		h_subArray[0] += h_subArray[i];
//	}
//
//	T e = h_addArray[0] - h_subArray[0];
//
//	return e;
//	}

void writeToFile(int numIter, int numBodies, double timeToMomentum, /*float energy1,*/ double timeToScale, double timeToAdvance, double timeToScale2, /*float energy2,*/ double total)
{
	file << numIter << "," << numBodies << "," << timeToScale << ',' << timeToScale2 << "," << timeToAdvance << "," << timeToMomentum << "," <</* energy1 << ',' << energy2 << ',' <<*/ total << "\n";
}

template <typename T>
void gpuLoops(int niters, int nbodies)
{
	
		type e1, e2;
		auto t1 = std::chrono::steady_clock::now();
		auto t2 = std::chrono::steady_clock::now();

		auto Tadv = std::chrono::steady_clock::now();
		auto Tadv2 = std::chrono::steady_clock::now();

		/*auto t1 = std::chrono::steady_clock::now();
		auto t2 = std::chrono::steady_clock::now();*/

		auto momentumStart = std::chrono::steady_clock::now();
		auto momentumEnd = std::chrono::steady_clock::now();

		//auto energyStart1 = std::chrono::steady_clock::now();
		//auto energyEnd1 = std::chrono::steady_clock::now();

		auto scaleStart1 = std::chrono::steady_clock::now();
		auto scaleEnd1 = std::chrono::steady_clock::now();

		auto advStart = std::chrono::steady_clock::now();
		auto advEnd = std::chrono::steady_clock::now();

		auto scaleStart2 = std::chrono::steady_clock::now();
		auto scaleEnd2 = std::chrono::steady_clock::now();

		planet<type> *bodies;
		if (nbodies == 5) {
			bodies = golden_bodies; // Check accuracy with 1000 solar system iterations
		}
		else {
			bodies = new planet<type>[nbodies];
			init_random_bodies(nbodies, bodies);
		}
		planet<type> *Gbodies;

		cudaMalloc(&Gbodies, nbodies*sizeof(planet<type>));
		cudaMemcpy(Gbodies, bodies, nbodies*sizeof(planet<type>), cudaMemcpyHostToDevice);

		int maxThreads;
		cudaDeviceGetAttribute(&maxThreads, cudaDevAttrMaxThreadsPerBlock, 0);
		if (nbodies < maxThreads)
			blockSize = nbodies;
		else
			blockSize = maxThreads;
		gridSize = (nbodies + blockSize) / blockSize;

		t1 = std::chrono::steady_clock::now();

		momentumStart = std::chrono::steady_clock::now();
		callOffSet(nbodies, Gbodies);
		momentumEnd = std::chrono::steady_clock::now();

		cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess)
		{
			std::cout << "Error in position kernal: " << cudaGetErrorString(error) << std::endl;
		}
		cudaThreadSynchronize();

		cudaMemcpy(bodies, Gbodies, nbodies*sizeof(planet<type>), cudaMemcpyDeviceToHost);
		e1 = energy(nbodies, bodies);
		cudaMemcpy(Gbodies, bodies, nbodies*sizeof(planet<type>), cudaMemcpyHostToDevice);

		//Scaling initial
		scaleStart1 = std::chrono::steady_clock::now();
		scale_bodies_GPU << <gridSize, blockSize >> >(nbodies, Gbodies, DT);
		scaleEnd1 = std::chrono::steady_clock::now();


		//Tadv = std::chrono::steady_clock::now();
		advStart = std::chrono::steady_clock::now();
		for (auto i = 0; i < niters; ++i)
		{
			//Calling advanced
			adv_Update_GPU << <gridSize, blockSize >> >(nbodies, Gbodies);
			cudaDeviceSynchronize();

			//	adv_Position_Update << <gridSize, blockSize >> >(nbodies, Gbodies);
		}
		advEnd = std::chrono::steady_clock::now();

		//Scaling again
		scaleStart2 = std::chrono::steady_clock::now();
		scale_bodies_GPU << <gridSize, blockSize >> >(nbodies, Gbodies, RECIP_DT);
		scaleEnd2 = std::chrono::steady_clock::now();

		cudaMemcpy(bodies, Gbodies, nbodies*sizeof(planet<type>), cudaMemcpyDeviceToHost);
		e2 = energy(nbodies, bodies);

		t2 = std::chrono::steady_clock::now();

		auto momDiff = momentumEnd - momentumStart;
		auto TimeMom = std::chrono::duration<double>(momDiff).count();

		//auto energyDiff = Tenergy2 - Tenergy;
		//auto  TimeEnergy = std::chrono::duration<double>(energyDiff).count();

		auto scaleDiff = scaleEnd1 - scaleStart1;
		auto TimeSc = std::chrono::duration<double>(scaleDiff).count();

		auto advDiff = advEnd - advStart;
		auto TimeAdv = std::chrono::duration<double>(advDiff).count();

		auto scale2Diff = scaleEnd2 - scaleStart2;
		auto TimeSc2 = std::chrono::duration<double>(scale2Diff).count();

		//auto energy2Diff = Tenergy2 - Tenergy;
		//auto  TimeEnergy2 = std::chrono::duration<double>(energy2Diff).count();

		auto diff = t2 - t1;
		auto TimeTotal = std::chrono::duration<double>(diff).count();

		writeToFile(niters, nbodies, TimeMom, /*energy1T*/ TimeSc, TimeAdv, TimeSc, /*energy2T*/ TimeTotal);
		std::cout << "part done \n";
}

int main(int argc, char ** argv)
{
	int niters = 1000, nbodies = 900;
	if (argc > 1) { niters = atoi(argv[1]); }
	if (argc > 2) { nbodies = atoi(argv[2]); }

	std::cout << "niters=" << niters << " nbodies=" << nbodies << '\n';

	outData = new type[gridSize];


	//if (!GPUTEST)
	//{
	//	t1 = std::chrono::steady_clock::now();
	//	offset_momentum(nbodies, bodies);
	//	e1 = energy(nbodies, bodies);
	//	scale_bodies(nbodies, bodies, DT);
	//	for (int i = 1; i <= niters; ++i)  {
	//		advance(nbodies, bodies);
	//	}
	//	scale_bodies(nbodies, bodies, RECIP_DT);

	//	e2 = energy(nbodies, bodies);
	//	t2 = std::chrono::steady_clock::now();
	//}
	
		file.open("Test.csv");
		file << "Iterations" << ',' << "Body Count" << ',' << "Time for Scale" << ',' << "Time for Scale 2" << ',' << "Time for Advance" << ',' << "Time for Momentum" << ',' <</* "Energy Before" << ',' << "Energy After" << ',' <<*/ "Total" << '\n';
		for (nbodies = 100; nbodies <= 1000; nbodies += 100)
		{
			for (niters = 100; niters <= 1000; niters += 100)
			{
	
				gpuLoops<type>(niters, nbodies);
			}
		}
		file.close();
		//Free up the memory
	/*	cudaFree(Gbodies);
	}*/

	/*auto diff = t2 - t1;
	auto diff2 = Tadv2 - Tadv;

	std::cout << std::setprecision(9);
	std::cout << e1 << '\n' << e2 << '\n';
	std::cout << std::chrono::duration<double>(diff).count() << " seconds.\n";
	std::cout <<"adv: " << std::chrono::duration<double>(diff2).count() << " seconds.\n";
	std::cout << "GridSize: " << gridSize << std::endl;*/
	std::cout << "DONE \n";
	std::cin.get();

	delete[]outData;
	//if (argc != 1) { delete[] bodies; }
	return 0;
}
