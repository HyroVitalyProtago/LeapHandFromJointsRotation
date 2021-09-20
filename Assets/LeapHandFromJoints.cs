using System;
using System.Collections;
using System.Collections.Generic;
using Leap;
using Leap.Unity;
using Leap.Unity.Attributes;
using Leap.Unity.Infix;
using UnityEngine;

public class LeapHandFromJoints : MonoBehaviour {

	private const float CYLINDER_MESH_RESOLUTION = 0.1f; //in centimeters, meshes within this resolution will be re-used
	[MinValue(3)]
	[SerializeField]
	private int _cylinderResolution = 12;
	
	[MinValue(0)]
	[SerializeField]
	private float _cylinderRadius = 0.006f;
	
	const int NUM_POINTS = 22;
	
	static readonly double[][] points = {
		new[] {
			26.89253807, 50.71062469, 26.30517006,
			-6.53940153, -28.20298195, -43.96928406,
			1.59295416, -18.24504852, -15.59684658,
			-9.10624218, 18.74389076, -4.08834696,
			-18.99320984, -28.7378273, 36.90289307,
			19.9671917, 6.81442356, -2.79738522,
			51.97451019, 43.90403748, 36.57212448, 28.60396194
		},
		new[] {
			35.12473297, 67.64217377, 79.15477753, 69.94815826, 60.34004593,
			58.72478104, 32.09606171, 11.13847065, 16.22161102, 24.88112068,
			22.31115341, -14.11880112, -31.75012016, -40.87027359, 15.9162302,
			-20.30714417, -39.99034882, -51.30535889, 10.07509327, -20.53441429,
			-36.16606903, -48.56869125
		},
		new[] {
			151.66641235, 123.82020569, 118.66936493, 147.56196594, 166.81466675,
			180.52867126, 170.5990448, 145.15351868, 124.28096008, 113.46463776,
			168.91151428, 165.07273865, 154.32685852, 144.07644653, 162.14595032,
			160.73678589, 153.32025146, 145.6038208, 151.8303833, 152.60569763,
			149.65892029, 145.08146667
		}
	};

	static readonly double[][][] angles = new[] {
		new[] {
			new[] {6.80335651e-09, -4.51514916e-01, -1.10714834e+00},
			new[] {-2.12284823e-01, -1.89584131e-02, 4.08588116e-03},
			new[] {1.12782320e-01, -3.36439092e-09, -3.76808121e-09},
			new[] {-2.43190810e-01, -1.47738551e-08, -3.13275023e-09}
		},

		new[] {
			new[] {1.51001870e-01, -1.66816955e-01, -2.52599157e-02},
			new[] {1.26234886e+00, 1.18541579e-03, 3.72051055e-03},
			new[] {1.10766007e+00, -9.07727204e-09, 1.62977681e-08},
			new[] {5.22654522e-01, -1.12897120e-08, 1.53981473e-08}
		},

		new[] {
			new[] {1.48955348e-01, -2.95250742e-02, 1.46350714e-01},
			new[] {5.97808271e-01, 1.37448961e-01, 9.30286183e-02},
			new[] {3.72391062e-01, -4.33029963e-08, -9.43238189e-09},
			new[] {2.40749485e-01, -2.57356516e-08, -1.51517855e-08}
		},

		new[] {
			new[] {1.50006118e-01, 1.21617516e-01, 2.20302218e-01},
			new[] {5.11299735e-01, 1.47075129e-01, 8.20371808e-02},
			new[] {3.06999060e-01, 2.59781033e-08, -1.21896506e-08},
			new[] {2.05311882e-01, 3.33055462e-08, -4.01401060e-09}
		},

		new[] {
			new[] {1.09821361e-01, 2.62327266e-01, 3.74564087e-01},
			new[] {4.06376774e-01, 2.33796594e-01, 9.93673235e-02},
			new[] {2.64315977e-01, 2.51456589e-08, 2.37648319e-08},
			new[] {1.84392151e-01, -5.94240413e-09, 5.29384378e-09}
		}
	};

	[SerializeField]
	private bool _castShadows = true;

	[SerializeField]
	private Material _material;
    
	[SerializeField]
	private Mesh _sphereMesh;
    
	[MinValue(0)]
	[SerializeField]
	private float _jointRadius = 0.008f;
    
	[MinValue(0)]
	[SerializeField]
	private float _palmRadius = 0.015f;
    
	private const int TOTAL_JOINT_COUNT = 4 * 5;
	private const int THUMB_BASE_INDEX = (int)Finger.FingerType.TYPE_THUMB * 4;
    
	LeapProvider provider;
	private Material _sphereMat;
	private Vector3[] _spherePositions = new Vector3[TOTAL_JOINT_COUNT];
    
	void Awake() {
		provider = GetComponent<LeapProvider>();
		if (_material != null) {
			_sphereMat = new Material(_material);
			_sphereMat.hideFlags = HideFlags.DontSaveInEditor;
		}

		//handGo = new GameObject();
	}

	/*
	[Range(0,5)]
	public int cse = 0;

	GameObject[] gos = new GameObject[5 * 5];
	GameObject handGo;
	void Update1() {
		// From positions
		for (int i = 0; i < NUM_POINTS; i++) {
			drawSphere(new Vector3((float)points[0][i]/1000f, (float)points[1][i]/1000f, (float)points[2][i]/1000f));
		}
	
		for (int finger = 0; finger < 5; finger++) {
			for (int bone = 0; bone < 4-1; bone++) {
				int i = 2 + finger*4 + bone;
				var p1 = new Vector3((float) points[0][i] / 1000f, (float) points[1][i] / 1000f,
					(float) points[2][i] / 1000f);
				int j = 2 + finger*4 + (bone + 1);
				var p2 = new Vector3((float) points[0][j] / 1000f, (float) points[1][j] / 1000f,
					(float) points[2][j] / 1000f);
				drawCylinder(p1,p2);
			}
		}
	
		var scale = .1f;
		var length = new float[] {0, .6f*scale, .3f*scale, .25f*scale, .2f*scale};
	
		var palmPos = new Vector3((float) points[0][0] / 1000f, (float) points[1][0] / 1000f,
			(float) points[2][0] / 1000f);
		var wristPos = new Vector3((float) points[0][1] / 1000f, (float) points[1][1] / 1000f,
			(float) points[2][1] / 1000f);
	
		handGo.transform.position = wristPos;
		
		// From rotations
		for (int finger = 0; finger < 5; finger++) {
			Vector3 pos = Vector3.zero;
			Quaternion rot = Quaternion.identity;
			GameObject go = handGo;
			for (int bone = 0; bone < 4; bone++) {
				var pitch = angles[finger][bone][0] * Mathf.Rad2Deg;
				var yaw = angles[finger][bone][1] * Mathf.Rad2Deg;
				var roll = angles[finger][bone][2] * Mathf.Rad2Deg;
	
				int id = finger * 5 + bone;
				if (gos[id] == null) gos[id] = new GameObject();
				var g = gos[id];
				if (go != null) g.transform.SetParent(go.transform);
	
				int i = 2 + finger*4 + bone-1;
				var p1 = new Vector3((float) points[0][i] / 1000f, (float) points[1][i] / 1000f,
					(float) points[2][i] / 1000f);
				int j = 2 + finger*4 + bone;
				var p2 = new Vector3((float) points[0][j] / 1000f, (float) points[1][j] / 1000f,
					(float) points[2][j] / 1000f);
	
	
				g.transform.localPosition = new Vector3(0, 0, finger == 0 && bone == 1 ? 0 : length[bone]);
				g.transform.localRotation = Quaternion.Euler((float) yaw, (float) pitch, (float) roll);
	
				
				var oldPos = pos;
				pos = g.transform.position;
				go = g;
				drawSphere(pos);
	
				if (bone > 0 && !(finger == 0 && bone == 1)) {
					drawCylinder(oldPos, pos);
				}
			}
			
			int id2 = finger * 5 + 4;
			if (gos[id2] == null) gos[id2] = new GameObject();
			var g2 = gos[id2];
			if (go != null) g2.transform.SetParent(go.transform);
	
			g2.transform.localPosition = new Vector3(0, 0, length[4]);
			drawCylinder(pos, g2.transform.position);
			drawSphere(g2.transform.position);
		}
	
	}

	void Update2() {
		var hands = provider.CurrentFrame.Hands;
		if (hands.Count <= 0) return;

		var hand = hands[0];
        
		Vector3 palmPosition = hand.PalmPosition.ToVector3();
		drawSphere(Vector3.zero, _palmRadius);
        
		// Update all joint spheres in the fingers
		foreach (var finger in hand.Fingers) {
			var prevPos = Vector3.zero;
			var prevRot = Quaternion.identity;
			
			GameObject prevJoint = null;
			for (int j = 0; j < 4; j++) {
				var bone = finger.Bone((Bone.BoneType) j);
				var length = bone.Length;
				
				// drawSphere(bone.PrevJoint.ToVector3());
				// drawSphere(bone.NextJoint.ToVector3());
				// if (length > 0) drawCylinder(bone.PrevJoint.ToVector3(), bone.NextJoint.ToVector3());
				
				int id = ((int)finger.Type) * 5 + j;
				if (gos[id] == null) gos[id] = new GameObject();
				var g = gos[id];
				
				 if (prevJoint != null) {
				 	g.transform.SetParent(prevJoint.transform);
				 }

				 g.transform.rotation = Quaternion.Inverse(hand.Rotation.ToQuaternion()) * bone.Rotation.ToQuaternion();
				 g.transform.position = bone.Center.ToVector3();

				 var rot = prevRot * g.transform.localRotation;
				 if (j == 0) prevPos = Quaternion.Inverse(hand.Rotation.ToQuaternion()) * (bone.PrevJoint.ToVector3() - hand.PalmPosition.ToVector3());
				 var nextPos = prevPos + rot * new Vector3(0, 0, length);
				 drawSphere(prevPos);
				 drawSphere(nextPos);
				 if (length > 0) drawCylinder(prevPos, nextPos);

				 prevPos = nextPos;
				 prevRot = rot;
				 
				// drawSphere(g.transform.position);
				// if (prevJoint != null) {
				// 	drawCylinder(prevJoint.transform.position, g.transform.position);
				// }
				
				prevJoint = g;
			}
		}
	}
*/
	void Update3() {
		var hands = provider.CurrentFrame.Hands;
		if (hands.Count <= 0) return;

		var hand = hands[0];

		var angles = GetAngles(hand);
		for (int finger = 0; finger < 5; finger++) { // finger
			var prevPos = Quaternion.Inverse(hand.Rotation.ToQuaternion()) * (hand.Fingers[finger].bones[0].PrevJoint.ToVector3() - hand.PalmPosition.ToVector3());//Vector3.zero;
			var prevRot = Quaternion.identity;
			var prevRotMat = Matrix.Identity;
			for (int bone = 0; bone < 4; bone++) { // bones
				var pitch = angles[finger, bone, 0] * Mathf.Rad2Deg;
				var yaw = angles[finger, bone, 1] * Mathf.Rad2Deg;
				var roll = angles[finger, bone, 2] * Mathf.Rad2Deg;

				//Quaternion rot = prevRot * Quaternion.Euler((float) pitch, (float) -yaw, (float) roll);
				//rot = prevRot * ToQuaternion(yaw, pitch, roll);

				Matrix rotMat = prevRotMat * GetRotFromAngles(pitch * Mathf.Deg2Rad, yaw * Mathf.Deg2Rad, roll * Mathf.Deg2Rad);
				Quaternion rot = prevRot * FromRotationMat(rotMat);
				
				//if (bone == 0) prevPos = Quaternion.Inverse(hand.Rotation.ToQuaternion()) * (hand.Fingers[finger].bones[bone].PrevJoint.ToVector3() - hand.PalmPosition.ToVector3());
				
				//var nextPos = prevPos + rot * new Vector3(0, 0, hand.Fingers[finger].bones[bone].Length);
				var nextPos = prevPos + rotMat.TransformDirection(new Vector(0, 0, hand.Fingers[finger].bones[bone].Length)).ToVector3();

				drawSphere(prevPos);
				drawSphere(nextPos);
				if (finger != 0 || bone != 0) drawCylinder(prevPos, nextPos);

				prevPos = nextPos;
				prevRot = rot;
				prevRotMat = rotMat;
			}
		}
	}
	
	Quaternion ToQuaternion(float yaw, float pitch, float roll) // yaw (Z), pitch (Y), roll (X)
	{
		yaw *= Mathf.Deg2Rad;
		pitch *= Mathf.Deg2Rad;
		roll *= Mathf.Deg2Rad;
		
		// Abbreviations for the various angular functions
		float cy = Mathf.Cos(yaw * 0.5f);
		float sy = Mathf.Sin(yaw * 0.5f);
		float cp = Mathf.Cos(pitch * 0.5f);
		float sp = Mathf.Sin(pitch * 0.5f);
		float cr = Mathf.Cos(roll * 0.5f);
		float sr = Mathf.Sin(roll * 0.5f);

		Quaternion q = new Quaternion(
			w: cr * cp * cy + sr * sp * sy,
			x: sr * cp * cy - cr * sp * sy,
			y: cr * sp * cy + sr * cp * sy,
			z: cr * cp * sy - sr * sp * cy);

		return q;
	}

	// https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
	Quaternion FromRotationMat(Matrix m) {
		var ma = m.ToArray3x3();
		var m00 = ma[0];
		var m01 = ma[1];
		var m02 = ma[2];
		var m10 = ma[3];
		var m11 = ma[4];
		var m12 = ma[5];
		var m20 = ma[6];
		var m21 = ma[7];
		var m22 = ma[8];
		
		Quaternion q = Quaternion.identity;
		float t;
		
		if (m22 < 0) {
			if (m00 > m11) {
				t = 1 + m00 - m11 - m22;
				q = new Quaternion(t, m01 + m10, m20 + m02, m12 - m21);
			} else {
				t = 1 - m00 + m11 - m22;
				q = new Quaternion(m01 + m10, t, m12 + m21, m20 - m02);
			}
		} else {
			if (m00 < -m11) {
				t = 1 - m00 - m11 + m22;
				q = new Quaternion(m20 + m02, m12 + m21, t, m01 - m10);
			} else {
				t = 1 + m00 + m11 + m22;
				q = new Quaternion(m12 - m21, m20 - m02, m01 - m10, t);
			}
		}

		return ScalarMultiply(q, .5f / Mathf.Sqrt(t));
	}
	
	static Quaternion ScalarMultiply(Quaternion quaternion, float scalar) {
		return new Quaternion(quaternion.x * scalar, quaternion.y * scalar, quaternion.z * scalar, quaternion.w * scalar);
	}
	
	Matrix GetRotationMat(LeapTransform lt) {
		return new Matrix(lt.xBasis, lt.yBasis, lt.zBasis);
	}

	Vector3 AnglesFromRot(Matrix mat) {
		/*
		Function from LearnOpenCV, Satya Mallick:
		https://www.learnopencv.com/rotation-matrix-to-euler-angles/
		https://github.com/spmallick/learnopencv/blob/master/RotationMatrixToEulerAngles/rotm2euler.py
		*/
		var sy = Mathf.Sqrt(mat.xBasis[0] * mat.xBasis[0] + mat.yBasis[0] * mat.yBasis[0]);

		if (sy < 1e-6) {
			return new Vector3(
				Mathf.Atan2(mat.zBasis[1], mat.zBasis[2]),
				Mathf.Atan2(-mat.zBasis[0], sy),
				Mathf.Atan2(mat.yBasis[0], mat.xBasis[0])
				);
		}
		return new Vector3(
			Mathf.Atan2(-mat.yBasis[2], mat.yBasis[1]),
			Mathf.Atan2(-mat.zBasis[0], sy),
			0
		);
	}

	/*
	 Function from LearnOpenCV, Satya Mallick:
	 https://www.learnopencv.com/rotation-matrix-to-euler-angles/
	 https://github.com/spmallick/learnopencv/blob/master/RotationMatrixToEulerAngles/rotm2euler.py
	 */
	Matrix GetRotFromAngles(float pitch, float yaw, float roll) {
		var x = new Matrix(
			1, 0,                0,
			0, Mathf.Cos(pitch), -Mathf.Sin(pitch),
			0, Mathf.Sin(pitch), Mathf.Cos(pitch)
		);
		
		var y = new Matrix(
			Mathf.Cos(yaw),  0, Mathf.Sin(yaw),
			0,               1, 0,
			-Mathf.Sin(yaw), 0, Mathf.Cos(yaw)
		);
		
		var z = new Matrix(
			Mathf.Cos(roll), -Mathf.Sin(roll), 0,
			Mathf.Sin(roll), Mathf.Cos(roll),  0,
			0,               0,                1
		);

		return z * (y * x);
	}
	
	float[,,] GetAngles(Hand hand, bool drawGizmos = false) {
		float[,,] ret = new float[5,4,3];
		for (int i = 0; i < 5; i++) {
			Finger finger = hand.Fingers[i];
			for (int j = 0; j < 4; j++) {
				Bone bone = finger.bones[j];

				// Rotation matrix from Basis vector
				var lastBoneMat = /*GetRotationMat(hand.Basis).RigidInverse() */ GetRotationMat(j == 0 ? hand.Basis : finger.bones[j - 1].Basis);
				var boneMat = /*GetRotationMat(hand.Basis).RigidInverse() */ GetRotationMat(bone.Basis);

				if (drawGizmos) {
					DrawAxis(bone.Center.ToVector3(), boneMat.xBasis.ToVector3(), boneMat.yBasis.ToVector3(), boneMat.zBasis.ToVector3());
				}

				// Rotation matrix between bones
				var rotMat = lastBoneMat.RigidInverse() * boneMat;
				
				// Euler angles from rotation matrix
				var anglesFromRot = AnglesFromRot(rotMat);

				ret[i, j, 0] = anglesFromRot.x;
				ret[i, j, 1] = anglesFromRot.y;
				ret[i, j, 2] = anglesFromRot.z;
			}
		}

		return ret;
	}

	Matrix Transpose(Matrix m) {
		var ma = m.ToArray3x3();
		return new Matrix(ma[0], ma[3], ma[6], ma[1], ma[4], ma[7], ma[2], ma[5], ma[8]);
	}

	void DrawAxis(Vector3 from, Vector3 x, Vector3 y, Vector3 z) {
		float length = .01f;
		
		Gizmos.color = Color.red;
		Gizmos.DrawLine(from, from + (x - from)*length);
					
		Gizmos.color = Color.green;
		Gizmos.DrawLine(from, from + (y - from)*length);
					
		Gizmos.color = Color.blue;
		Gizmos.DrawLine(from, from + (z - from)*length);
	}
	
	void Update() => Update3();

	void OnDrawGizmos() {
		if (provider == null || provider.CurrentFrame == null) return;
		var hands = provider.CurrentFrame.Hands;
		if (hands.Count <= 0) return;

		var hand = hands[0];

		var angles = GetAngles(hand, true);
	}

	private int getFingerJointIndex(int fingerIndex, int jointIndex) {
		return fingerIndex * 4 + jointIndex;
	}
    
	private void drawSphere(Vector3 position) {
		drawSphere(position, _jointRadius);
	}
    
	private void drawSphere(Vector3 position, float radius) {
		//multiply radius by 2 because the default unity sphere has a radius of 0.5 meters at scale 1.
		Graphics.DrawMesh(_sphereMesh, 
			Matrix4x4.TRS(position, 
				Quaternion.identity, 
				Vector3.one * radius * 2.0f * transform.lossyScale.x), 
			_sphereMat, 0, 
			null, 0, null, true);
	}
	
	 private void drawCylinder(Vector3 a, Vector3 b) {
      float length = (a - b).magnitude;

      Graphics.DrawMesh(getCylinderMesh(length),
                        Matrix4x4.TRS(a, 
                                      Quaternion.LookRotation(b - a), 
                                      new Vector3(transform.lossyScale.x, transform.lossyScale.x, 1)),
                        _material,
                        gameObject.layer, 
                        null, 0, null, _castShadows);
    }

    private void drawCylinder(int a, int b) {
      drawCylinder(_spherePositions[a], _spherePositions[b]);
    }

    private void drawCylinder(Vector3 a, int b) {
      drawCylinder(a, _spherePositions[b]);
    }

    private Dictionary<int, Mesh> _meshMap = new Dictionary<int, Mesh>();
    private Mesh getCylinderMesh(float length) {
      int lengthKey = Mathf.RoundToInt(length * 100 / CYLINDER_MESH_RESOLUTION);

      Mesh mesh;
      if (_meshMap.TryGetValue(lengthKey, out mesh)) {
        return mesh;
      }

      mesh = new Mesh();
      mesh.name = "GeneratedCylinder";
      mesh.hideFlags = HideFlags.DontSave;

      List<Vector3> verts = new List<Vector3>();
      List<Color> colors = new List<Color>();
      List<int> tris = new List<int>();

      Vector3 p0 = Vector3.zero;
      Vector3 p1 = Vector3.forward * length;
      for (int i = 0; i < _cylinderResolution; i++) {
        float angle = (Mathf.PI * 2.0f * i) / _cylinderResolution;
        float dx = _cylinderRadius * Mathf.Cos(angle);
        float dy = _cylinderRadius * Mathf.Sin(angle);

        Vector3 spoke = new Vector3(dx, dy, 0);

        verts.Add(p0 + spoke);
        verts.Add(p1 + spoke);

        colors.Add(Color.white);
        colors.Add(Color.white);

        int triStart = verts.Count;
        int triCap = _cylinderResolution * 2;

        tris.Add((triStart + 0) % triCap);
        tris.Add((triStart + 2) % triCap);
        tris.Add((triStart + 1) % triCap);

        tris.Add((triStart + 2) % triCap);
        tris.Add((triStart + 3) % triCap);
        tris.Add((triStart + 1) % triCap);
      }

      mesh.SetVertices(verts);
      mesh.SetIndices(tris.ToArray(), MeshTopology.Triangles, 0);
      mesh.RecalculateBounds();
      mesh.RecalculateNormals();
      mesh.UploadMeshData(true);

      _meshMap[lengthKey] = mesh;

      return mesh;
    }
}
