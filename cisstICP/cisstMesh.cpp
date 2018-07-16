// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Ayushi Sinha, Russell Taylor, 
//	  Johns Hopkins University. 
//	  All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// ****************************************************************************

#include "cisstMesh.h"
#include "utilities.h"

#include <cisstNumerical/nmrLSSolver.h>

#include <assert.h>
#undef NDEBUG       // enable assert in release mode

#include <fstream>

// compares vertex for storing in std::Map
//  NOTE: this routine is used for loading a mesh from an STL file;
//        it is used to detect multiple copies of the same vertex
//        in order to prevent storing the same vertex coordinate to
//        multiple locations in the vertices array
struct VertexCompare {
  bool operator() (const vct3 &k1, const vct3 &k2) const
  { // Return true if k1 goes before k2 in the strict weak 
    //  ordering for the map object;
    //  i.e. if k1 is not strictly less than k2, return false
    //       otherwise, return true
    //
    // Note: https://ece.uwaterloo.ca/~dwharder/aads/Relations/Weak_ordering/

    // order by x,y,z in that order
    if (k1.X() < k2.X())
      return true;
    else if (k1.X() == k2.X()) {
      if (k1.Y() < k2.Y())
        return true;
      else if (k1.Y() == k2.Y()) {
        if (k1.Z() < k2.Z())
          return true;
      }
    }
    return false;
  };
};


void cisstMesh::ResetMesh()
{
  vertices.SetSize(0);
  faces.SetSize(0);
  faceNormals.SetSize(0);
  vertexNormals.SetSize(0);
  faceNeighbors.SetSize(0);
  TriangleCov.SetSize(0);
  TriangleCovEig.SetSize(0);
}

void cisstMesh::ResetModel()
{
	meanShape.SetSize(0);
	mode.SetSize(0);
	modeWeight.SetSize(0);
	wi.SetSize(0);
	Si.SetSize(0);
	//estVertices.SetSize(0);
}

void cisstMesh::InitializeNoiseModel()
{
  TriangleCov.SetSize(faces.size());
  TriangleCovEig.SetSize(faces.size());

  TriangleCov.SetAll(vct3x3(0.0));
  TriangleCovEig.SetAll(vct3(0.0));
}

void cisstMesh::InitializeNoiseModel(
  double noiseInPlaneVar,
  double noisePerpPlaneVar)
{
  vct3x3 M, M0;
  vctRot3 R;
  vct3 z(0.0, 0.0, 1.0);
  vct3 norm;

  if (faceNormals.size() != faces.size())
  {
    std::cout << "ERROR: must initialize face normals in order to compute mesh noise model" << std::endl;
    TriangleCov.SetSize(0);
    TriangleCovEig.SetSize(0);
    assert(0);
  }

  TriangleCov.SetSize(faces.size());
  TriangleCovEig.SetSize(faces.size());

  // set covariance eigenvalues (in order of decreasing magnitude)
  if (noiseInPlaneVar >= noisePerpPlaneVar)
  {
    TriangleCovEig.SetAll(vct3(noiseInPlaneVar, noiseInPlaneVar, noisePerpPlaneVar));
  }
  else
  {
    TriangleCovEig.SetAll(vct3(noisePerpPlaneVar, noiseInPlaneVar, noiseInPlaneVar));
  }

  // compute covariance matrices
  for (unsigned int i = 0; i < faces.size(); i++)
  {
    TriangleCov[i] = ComputePointCovariance(faceNormals[i], noisePerpPlaneVar, noiseInPlaneVar);
  }
}

void cisstMesh::SaveTriangleCovariances(std::string &filePath)
{
  std::cout << "Saving mesh covariances to file: " << filePath << std::endl;
  std::ofstream fs(filePath.c_str());
  if (!fs.is_open())
  {
    std::cout << "ERROR: failed to open file for saving cov: " << filePath << std::endl;
    assert(0);
  }
  unsigned int numCov = this->TriangleCov.size();
  //fs << "NUMCOV " << numCov << "\n";
  for (unsigned int i = 0; i < numCov; i++)
  {
    fs << this->TriangleCov.at(i).Row(0) << " "
      << this->TriangleCov.at(i).Row(1) << " "
      << this->TriangleCov.at(i).Row(2) << "\n";
  }
}

void cisstMesh::LoadPLY(const std::string &input_file) {
  ply_obj.read_ply(input_file,
    vertices, faces, faceNormals, faceNeighbors, vertexNormals);
  }

void cisstMesh::SavePLY(const std::string &output_file) {
  ply_obj.write_ply(output_file,
    vertices, faces, faceNormals, faceNeighbors, vertexNormals);
}

int cisstMesh::LoadMesh(
  const vctDynamicVector<vct3> *vertices,
  const vctDynamicVector<vctInt3> *faces,
  const vctDynamicVector<vct3> *face_normals,
  const vctDynamicVector<vctInt3> *face_neighbors,
  const vctDynamicVector<vct3> *vertex_normals
  )
{
  ResetMesh();

  if (!vertices || !faces) {
    std::cout << "ERROR: vertices and faces must not be null" << std::endl;
    assert(0);
  }
  if (vertices->size() < 1 || faces->size() < 1)
{
    std::cout << "ERROR: vertices and faces must not be empty" << std::endl;
    assert(0);
}
  this->vertices = *vertices;
  this->faces = *faces;

  if (face_normals) {
    if (face_normals->size() != this->faces.size()) {
      std::cout << "ERROR: number of face normals does not equal number of faces" << std::endl;
      assert(0);
  }
    this->faceNormals = *face_normals;
}

  if (face_neighbors) {
    if (face_neighbors->size() != this->faces.size()) {
      std::cout << "ERROR: number of face neighbors does not equal number of faces" << std::endl;
      assert(0);
  }
    this->faceNeighbors = *face_neighbors;
  }

  if (vertex_normals) {
    if (vertex_normals->size() != this->vertices.size()) {
      std::cout << "ERROR: number of face neighbors does not equal number of faces" << std::endl;
      assert(0);
  }
    this->vertexNormals = *vertex_normals;
    }

  InitializeNoiseModel();

  return 0;
}

int cisstMesh::LoadModelFile(const std::string &modelFilePath, int numModes)
{
	int rv;

	ResetModel();

	rv = AddModelFile(modelFilePath, numModes);

	return rv;
}

int cisstMesh::AddModelFile(const std::string &modelFilePath, int modes)
{
	//Load model from file having format:

	//file_location Nvertices=nvertices Nmodes=nmodes
	//Mode 0 : Mean Vertex Values
	//mx my mz
	//	...
	//mx my mz

	//Mode 1 : Vertex Displacements modeweight
	//mx my mz
	//	...
	//mx my mz

	//...

	//Mode nmodes : Vertex Displacements modeweight
	//mx my mz
	//	...
	//mx my mz
	float f1, f2, f3;
	unsigned int itemsRead;
	std::string line;

	unsigned int mOffset = meanShape.size();
	unsigned int wOffset = modeWeight.size();

	// open file
	std::ifstream modelFile;
	modelFile.open(modelFilePath.c_str());
	if (!modelFile.is_open())
	{
		std::cout << "ERROR: failed to open model file" << std::endl;
		return -1;
	}

	// read modes
	char fileLocation[100];
	unsigned int numVertices, numModes, modeNum;
	float modeWt;
	std::getline(modelFile, line);
	itemsRead = std::sscanf(line.c_str(), "%s Nvertices= %u Nmodes= %u", fileLocation, &numVertices, &numModes);
	if (itemsRead != 3)
	{
		std::cout << "ERROR: expected header at line: " << line << std::endl;
		return -1;
	}
	std::cout << " out of " << numModes ;
	if (vertices.size() != numVertices)
	{
		std::cout << "ERROR: model data does not match mesh data - number of vertices are different." << std::endl;
		return 0;
	}

	unsigned int modeCount = 0;
	vctDynamicVector<float> W;
	vctDynamicVector<float> V;
	vctDynamicVector<float> M;
	vctDynamicVector<float> E;
	W.SetSize(numVertices * 3);
	V.SetSize(numVertices * 3);
	M.SetSize(numVertices * 3);
	E.SetSize(numVertices * 3);

	vct3 m, w;
	modeWeight.resize(wOffset + modes - 1);
	meanShape.resize(mOffset + numVertices);
	mode.resize(wOffset + modes - 1);
	wi.resize(wOffset + modes - 1);
	Si.resize(mOffset + modes - 1);
	//estVertices.resize(mOffset + numVertices);

	while (modelFile.good() && modeCount < modes)
	{
		unsigned int vertCount = 0;
		std::getline(modelFile, line);
		if (modeCount < 1)
		{
			itemsRead = std::sscanf(line.c_str(), "Mode %u :Mean Vertex Values", &modeNum);
			if (itemsRead != 1)
			{
				std::cout << "ERROR: expected header at line: " << line << std::endl;
				return -1;
			}
			//std::cout << "\nMode : " << modeNum << "; Mean shape" << std::endl;
		}
		else
		{
			mode[wOffset + modeCount - 1].resize(mOffset + numVertices);
			wi[wOffset + modeCount - 1].resize(mOffset + numVertices);
			itemsRead = std::sscanf(line.c_str(), "Mode %u :Vertex Displacements %f", &modeNum, &modeWt);
			modeWeight[modeCount - 1] = modeWt;
			if (itemsRead != 2)
			{
				std::cout << "ERROR: expected header at line: " << line << std::endl;
				return -1;
			}
			//std::cout << "Mode : " << modeNum << "; Mode weight : " << modeWeight[modeCount - 1] << std::endl;
		}

		while (vertCount < numVertices)
		{
			std::getline(modelFile, line);
			itemsRead = std::sscanf(line.c_str(), "%f , %f , %f", &f1, &f2, &f3);
			if (itemsRead != 3)
			{
				std::cout << "ERROR: expected a point value at line: " << line << std::endl;
				return -1;
			}
			m[0] = f1;
			m[1] = f2;
			m[2] = f3;
			if (modeCount < 1)
				meanShape.at(mOffset + vertCount).Assign(m);
			else
			{
				mode[wOffset + modeCount - 1].at(mOffset + vertCount).Assign(m);

				// wi = sqrt(lambda_i)*mi
				w = m * sqrt(modeWeight[modeCount - 1]);
				wi[wOffset + modeCount - 1].at(mOffset + vertCount).Assign(w); 
			}
			vertCount++;
		}

		if (modelFile.bad() || modelFile.fail() || vertCount != numVertices)
		{
			std::cout << "ERROR: read points from model file failed; last line read: " << line << std::endl;
			return -1;
		}

		// Si = 0 (initialization of Si)
		if (modeCount > 0)
			// Set Si to 0
			Si[wOffset + modeCount - 1] = 0.0; 
		modeCount++;
	}

	// replace patient mesh with model estimate of the mesh
	this->vertices = meanShape; 

	if (modelFile.bad() || modelFile.fail() || modeCount != modes)
	{
		std::cout << "ERROR: read points from model file failed; last line read: " << line << std::endl;
		return -1;
	}

	std::cout << std::endl;
	return 1;
}