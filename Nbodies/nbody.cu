/*
   Running without arguments is equivalent to 1000 iterations with the
   5 celestial objects declared in the golden_bodies array.

   $ nbody.exe 1000 5

   The output of this shows the energy before and after the simulation,
   and should be:

   -0.169075164
   -0.169087605
*/

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#define GPUTEST 1


using type = double;

const type pi{3.141592653589793};
const type solar_mass{4 * pi * pi};
const type days_per_year{365.24};

template <typename T>
struct planet {
  T x, y, z;
  T vx, vy, vz;
  T mass;
};

//velocity update for the kernals
template <typename T>
__global__ void adv_Velocity_Update( int nbodies, planet<T> *bodies)
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
	}
}

//position update for the kernals
template <typename T>
__global__ void adv_Position_Update( int nbodies, planet<T> *bodies)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (i < nbodies)
	{
		planet<T> &b = bodies[i];
		b.x += b.vx;
		b.y += b.vy;
		b.z += b.vz;
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

template <typename T>
__global__ void energy_GPU(int nbodies, planet<T> *bodies)
{
	T e = 0.0;
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = i + 1;

	if (i < nbodies)
	{
		planet<T> &b = bodies[i];
		e += 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz);
		if (j < nbodies)
		{
			planet<T> &b2 = bodies[j];
			T dx = b.x - b2.x;
			T dy = b.y - b2.y;
			T dz = b.z - b2.z;
			T distance = sqrt(dx * dx + dy * dy + dz * dz);
			e -= (b.mass * b2.mass) / distance;
		}
	}
	//return e; need to return this type
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
      T inv_distance = 1.0/sqrt(dx * dx + dy * dy + dz * dz);
      T mag = inv_distance * inv_distance * inv_distance;
      b.vx  -= dx * b2.mass * mag;
      b.vy  -= dy * b2.mass * mag;
      b.vz  -= dz * b2.mass * mag;
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
  bodies[0].vx = - px / solar_mass;
  bodies[0].vy = - py / solar_mass;
  bodies[0].vz = - pz / solar_mass;
}

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

const type DT{1e-2};
const type RECIP_DT{1.0/DT};

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
    bodies[i].vx   *= scale;
    bodies[i].vy   *= scale;
    bodies[i].vz   *= scale;
  }
}

template <typename T>
void init_random_bodies(int nbodies, planet<T> *bodies)
{
  for (int i = 0; i < nbodies; ++i) {
    bodies[i].x    =  (T)rand()/RAND_MAX;
    bodies[i].y    =  (T)rand()/RAND_MAX;
    bodies[i].z    =  (T)rand()/RAND_MAX;
    bodies[i].vx   =  (T)rand()/RAND_MAX;
    bodies[i].vy   =  (T)rand()/RAND_MAX;
    bodies[i].vz   =  (T)rand()/RAND_MAX;
    bodies[i].mass =  (T)rand()/RAND_MAX;
  }
}

template <typename T>
void kernalUpdate(int nbodies, planet<T> *bodies)
{
	planet<T> *Gbodies;

	//Copy data from CPU
	cudaMalloc(&Gbodies, nbodies*sizeof(planet<T>));
	cudaMemcpy(Gbodies, bodies, nbodies, cudaMemcpyHostToDevice);

	//Scaling
	cudaMemcpy(Gbodies, bodies, nbodies*sizeof(planet<type>), cudaMemcpyHostToDevice);
	scale_bodies_GPU << <nbodies, ceil(nbodies / 2) >> >(nbodies, Gbodies, DT);
	cudaMemcpy(bodies, Gbodies, nbodies*sizeof(planet<type>), cudaMemcpyDeviceToHost);

	//velocity
	cudaMemcpy(Gbodies, bodies, nbodies*sizeof(planet<type>), cudaMemcpyHostToDevice);
	adv_Velocity_Update<<<1, nbodies>>>(nbodies, Gbodies);
	cudaMemcpy(bodies, Gbodies, nbodies*sizeof(planet<type>), cudaMemcpyDeviceToHost);

	//position
	cudaMemcpy(Gbodies, bodies, nbodies*sizeof(planet<type>), cudaMemcpyHostToDevice);
	adv_Position_Update << <nbodies, ceil(nbodies/2)>> >(nbodies, Gbodies);
	cudaMemcpy(bodies, Gbodies, nbodies*sizeof(planet<type>), cudaMemcpyDeviceToHost);

	//Scaling
	cudaMemcpy(Gbodies, bodies, nbodies*sizeof(planet<type>), cudaMemcpyHostToDevice);
	scale_bodies_GPU << <nbodies, ceil(nbodies / 2) >> >(nbodies, Gbodies, RECIP_DT);
	cudaMemcpy(bodies, Gbodies, nbodies*sizeof(planet<type>), cudaMemcpyDeviceToHost);

	//copy data back to CPU
	//cudaMemcpy(bodies, Gbodies, nbodies, cudaMemcpyDeviceToHost);

	//Free up the memory
	cudaFree(Gbodies);
}

int main(int argc, char ** argv)
{
  int niters = 1000, nbodies = 5;
  if (argc > 1) { niters  = atoi(argv[1]); }
  if (argc > 2) { nbodies = atoi(argv[2]); }

  std::cout << "niters=" << niters << " nbodies=" << nbodies << '\n';

  planet<type> *bodies;
  if (argc == 1) { 
    bodies = golden_bodies; // Check accuracy with 1000 solar system iterations
  } else {
    bodies = new planet<type>[nbodies];
    init_random_bodies(nbodies, bodies);
  }

  auto t1 = std::chrono::steady_clock::now();
  offset_momentum(nbodies, bodies);
  type e1 = energy(nbodies, bodies);
 
  if (GPUTEST)
  {
	  for (auto i = 0; i < niters; ++i)
		  kernalUpdate(nbodies, bodies);
  }
	else
	{
		  scale_bodies(nbodies, bodies, DT);
		  for (int i = 1; i <= niters; ++i)  {
		    advance(nbodies, bodies);
		  }
		  scale_bodies(nbodies, bodies, RECIP_DT);
	}
	 
  type e2 = energy(nbodies, bodies);
  auto t2 = std::chrono::steady_clock::now();
  auto diff = t2 - t1;

  std::cout << std::setprecision(9);
  std::cout << e1 << '\n' << e2 << '\n';
  std::cout << std::fixed << std::setprecision(3);
  std::cout << std::chrono::duration<double>(diff).count() << " seconds.\n";

  if (argc != 1) { delete [] bodies; }
  return 0;
}
