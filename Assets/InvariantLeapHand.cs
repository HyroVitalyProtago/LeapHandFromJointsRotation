using System;
using System.Collections;
using System.Collections.Generic;
using Leap;
using Leap.Unity;
using Leap.Unity.Attributes;
using Leap.Unity.Infix;
using UnityEngine;

public class InvariantLeapHand : MonoBehaviour {

    
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
    }

    void Update() {
        var hands = provider.CurrentFrame.Hands;
        if (hands.Count <= 0) return;

        var hand = hands[0];
        
        Vector3 palmPosition = hand.PalmPosition.ToVector3();
        drawSphere(Vector3.zero, _palmRadius);
        
        // Update all joint spheres in the fingers
        foreach (var finger in hand.Fingers) {
            for (int j = 0; j < 4; j++) {
                int key = getFingerJointIndex((int)finger.Type, j);

                Vector3 position = finger.Bone((Bone.BoneType)j).NextJoint.ToVector3();
                _spherePositions[key] = position;

                drawSphere((position - palmPosition).RotatedBy(Quaternion.Inverse(hand.GetPalmPose().rotation)));
            }
        }
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
}
