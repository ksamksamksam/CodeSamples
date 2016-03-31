using UnityEngine;
using System.Collections;

public class BehaviorScript : MonoBehaviour {


    public GameObject target;
    public GameObject floor;
    public Transform LeaderIndicator;
    GameObject indicator;
    
    TargetScript targetScript;
    GUIScript guiScript;
    
    float g_fMaxSpeed = 3.0f;
    float g_fMaxAccel = 10.0f;
    float g_fKNeighborhood = 100.0f;
    float g_fOriKv = 1.0f;  // Orientation
    float g_fOriKp = 1.0f; 
    float g_fVelKv = 1.0f;  // Velocity 
    
    float g_fAgentRadius = 1.0f;
    float g_fKArrival = 1.0f;
    float g_fKDeparture = 1.0f;
    float g_fKAlignment = 1.0f;
    float g_fKNoise = 1.0f;
    float g_fKWander = 1.0f;
    float g_fKSeparation = 1.0f;
    float g_fKCohesion = 1.0f;
    float g_fKAvoid = 1.0f;
    float g_fTAvoid = 1.0f;
    Vector3 m_vWander = Vector3.zero;
    
    Vector3 force;
    Vector3 torque;
    Vector3 targetPos;
    
    string behavior;

    // Use this for initialization
    void Start () {
        behavior = "seek";
        targetScript = (TargetScript)target.GetComponent( typeof(TargetScript) );   // contains agents
        guiScript = (GUIScript)floor.GetComponent( typeof(GUIScript) ); 
    }
    
    
    // Given an actor, update its movement through forces and torques.
    // This function obtains a desired velocity.  With the desired velocity
    // we can then calculate the corresponding torque and force that will
    // achieve it over time
    void FixedUpdate () {
        targetPos = target.transform.position;
        Vector3 dvel = CalculateDesiredVelocity(behavior);
        
        //if the desired velocity exceeds the max speed, change it to
        //a vector with its current direction at the max speed
        if (dvel.magnitude > g_fMaxSpeed)
        {
              dvel.Normalize();
              dvel *= g_fMaxSpeed;
        }
    
        force = CalculateForce(dvel);
        torque = CalculateTorque(dvel);
   
        //add the force and torque to the actor motion
        rigidbody.AddForce (force);
        rigidbody.AddTorque(torque);

    }



    // Given an actor and desired velocity, calculate a corresponding force
    Vector3 CalculateForce(Vector3 dvel)
    {
        //f = m[K0*(Vd-V)]
        Vector3 force = rigidbody.mass*(dvel-rigidbody.velocity);
        return force;
    }


    // Given an actor and desired velocity, calculate a corresponding torque
    Vector3 CalculateTorque(Vector3 dvel)
    {    
        // 1. Get current rotation matrix
         Matrix4x4 currRot = MatrixFromQuaternion(rigidbody.rotation);  
    
        // 2. Construct desired rotation matrix 
        // (This corresponds to facing in the direction of our desired velocity)
        Quaternion newQuat = new Quaternion();
        newQuat.SetFromToRotation(Vector3.forward, dvel);
        Matrix4x4 newRot = MatrixFromQuaternion(newQuat);
    
        // 3. Compute the change in rotation matrix that will
        // rotate the actor towards our desired rotation
        Matrix4x4 finRot = newRot*currRot.transpose;
        
        // 4. Construct quat for axisangle
        Quaternion finQuat = QuaternionFromMatrix(finRot);
        float an = 0.0f;
        Vector3 theta = Vector3.zero;
        finQuat.ToAngleAxis(out an, out theta);
    
        // find torque
        //
        Vector3 t = Vector3.Cross(rigidbody.angularVelocity, rigidbody.inertiaTensorRotation*rigidbody.angularVelocity)+rigidbody.inertiaTensorRotation*(g_fOriKp*an*theta-g_fOriKv*rigidbody.angularVelocity);
        
        return new Vector3(0, t.y, 0);
    }

    //set behavior algorithm to use
    public void setBehavior(string beh)
    {
        behavior = beh;
    }
    
    //Given a 4x4 matrix, calculate the quaternion
    public  Quaternion QuaternionFromMatrix(Matrix4x4 m) {
        
        Quaternion q = new Quaternion();
        q.w = Mathf.Sqrt( Mathf.Max( 0, 1 + m[0,0] + m[1,1] + m[2,2] ) ) / 2; 
        q.x = Mathf.Sqrt( Mathf.Max( 0, 1 + m[0,0] - m[1,1] - m[2,2] ) ) / 2; 
        q.y = Mathf.Sqrt( Mathf.Max( 0, 1 - m[0,0] + m[1,1] - m[2,2] ) ) / 2; 
        q.z = Mathf.Sqrt( Mathf.Max( 0, 1 - m[0,0] - m[1,1] + m[2,2] ) ) / 2; 
        q.x *= Mathf.Sign( q.x * ( m[2,1] - m[1,2] ) );
        q.y *= Mathf.Sign( q.y * ( m[0,2] - m[2,0] ) );
        q.z *= Mathf.Sign( q.z * ( m[1,0] - m[0,1] ) );
        return q;
    }
    
    //Given a quaternion, calculate a 4x4 matrix
    public Matrix4x4 MatrixFromQuaternion(Quaternion q)
    {
        
        float term00 = 1- 2 * (q.y * q.y + q.z * q.z);
        float term01 = 2 * (q.x * q.y - q.w * q.z);
        float term02 = 2 * (q.x * q.z + q.w * q.y);
    
    
        float term10 = 2 * (q.x * q.y + q.w * q.z);
        float term11 = 1 - 2 * (q.x * q.z + q.z * q.z);
        float term12 = 2 * (q.y * q.z - q.w * q.z);
    
    
        float term20 = 2 * (q.x * q.z - q.w * q.y);
        float term21 = 2 * (q.y * q.z + q.w * q.x);
        float term22 = 1 - 2 * (q.x * q.x + q.y * q.y);
            
        Vector4 column1  = new Vector4 (term00, term10, term20,0);
        Vector4 column2  = new Vector4 (term01, term11, term21,0);
        Vector4 column3  = new Vector4 (term02, term12, term22,0);
        
        
        Matrix4x4 M = new Matrix4x4();
        M.SetColumn(0, column1);
        M.SetColumn(1, column2);
        M.SetColumn(2, column3);
        M.SetColumn(3, new Vector4(0, 0, 0, 1));   
         
        return  M;
    }
    
    //Calculate a desired velocity based on the given behavior algorithm
    Vector3 CalculateDesiredVelocity(string behavior)
    {
    
        if(behavior == "leader")
        {
            LeaderIndicator.position = new Vector3 (targetScript.agents[0].transform.position.x, targetScript.agents[0].transform.position.y + 4.0f, targetScript.agents[0].transform.position.z);  // Place the leader indicator slightly above the leading actor
            return Leader(targetPos);
        }
        else
            LeaderIndicator.position = new Vector3(-2,-2,-2); //If not using Leader behavior, place the leader indicator below the floor
            
        if(behavior == "alignment")
            return Alignment();
        
        if(behavior == "arrival")
            return Arrival(targetPos);
            
        if(behavior == "cohesion")
            return Cohesion();
                
        if(behavior == "departure")
            return Departure(targetPos);
            
        if(behavior == "flee")
            return Flee(targetPos);
        
        if(behavior == "flocking")
            return Flocking();
        
        if(behavior == "seek")
            return Seek(targetPos);
            
        if(behavior == "separation")
            return Separation();    
            
        else
            return Wander(); //If no behavior is specified, have the actors wander the board aimlessly
            
     
    }
    
    Vector3 Alignment()
    {   
        //Have all the actors generally move in the same direction as their close neighbors.
        //Each actor moves with the average velocity of its neighbors.
        
        Vector3 sumVel = Vector3.zero;
        Vector3 dist = Vector3.zero;
        float count = 0;
        
        //Add up the velocities of every other actor within the distance defined by g_fKNeighborhood
        for (int i = 0; i < targetScript.agentNumber; i++)
        {
            dist = (rigidbody.position-targetScript.agents[i].rigidbody.position);
            if (dist.magnitude < g_fKNeighborhood && rigidbody != targetScript.agents[i].rigidbody)
            {
                count+=1;
                sumVel += targetScript.agents[i].rigidbody.velocity;
            }
        }
    
        //Return the average of these velocities, weighted by g_fKAlignment
        return g_fKAlignment*(sumVel/count);
    
    }
    
    Vector3 Arrival(Vector3 targetP)
    {
        //same behavior as seek except that the actor slows down once he is     
        //within range of the target
        
        //e = targetPosition - currentPosition
        Vector3 e = targetP - rigidbody.position;
        
        //check that actor is not already at target
        float dist = e.magnitude;
        //if he is, return the zero vector.  actor need not move.
        if (dist == 0)
            return Vector3.zero;
        
        //velocity damping factor, default set to 1, or no damping
        float damp = 1;
        //check if actor is within radius
        if (dist < g_fAgentRadius)
        {
            //if he is, damping factor is distance/radius
            damp = dist/g_fAgentRadius;
        }
        
        //Vd = K_arrival * e
        e *= g_fKArrival*damp;
        
        return e;
    }
    
    Vector3 Cohesion()
    {
        //Have the actors clump together as groups.
        //Each actor moves towards the centerpoint of its neighbors.
        
        Vector3 sumPos = Vector3.zero;
        float count = 0;
        
        //Add up the positions of every other actor within the distance defined by g_fKNeighborhood
        for (int i = 0; i < targetScript.agentNumber; i++)
        {
            Vector3 dist = (rigidbody.position-targetScript.agents[i].rigidbody.position);
            if (dist.magnitude < g_fKNeighborhood && rigidbody != targetScript.agents[i].rigidbody)
            {
                count+=1;
                sumPos += targetScript.agents[i].rigidbody.position;
            }
        }
        
        //Calculate the average position of the neighbors 
        Vector3 com = sumPos/count;
        
        //Use seeking behavior with the target being this centerpoint, weighted by g_fKCohesion
        return Seek(g_fKCohesion*(com-rigidbody.position) + rigidbody.position);
    }
    
    Vector3 Departure(Vector3 targetP)
    {
        //same behavior as flee except that actor only runs away 
        //if he is within radius and he runs slower as he gets farther
        
        //e = -(targetPosition - currentPosition)
        Vector3 e = -(targetP - rigidbody.position);
        
        float dist = e.magnitude;
        //check if actor is outside radius
        if (dist > g_fAgentRadius)
        {
            //if he is, he doesn't have to change velocity
            return rigidbody.velocity;
        }
        
        //velocity speedup factor, set to radius/distance
        float speedup = g_fAgentRadius/dist;
        
        //Vd = K_departure * e * speedup
        e *= g_fKDeparture*speedup;
        
        return e;
    }
    
    Vector3 Flee(Vector3 targetP)
    {
        //We want the actor to go as fast as he can away from the target
        //Vd = maxspeed * normalized_e
        
        //e = -(targetPosition - currentPosition)
        Vector3 e = -(targetP-rigidbody.position);
        
        //e = ||e|| * maxspeed
        e = e.normalized * g_fMaxSpeed;
        
        return e;
    }
    
    Vector3 Flocking()
    {
        //Model flocking as a weighted combination of separation,
        //cohesion, and alignment behavior.
        
        Vector3 Vflock = 8*Separation()+8*Cohesion()+10*Alignment();
      
        //Have the actor move at max speed in the calculated direction
        return g_fMaxSpeed*Vflock.normalized;
    }
    
    Vector3 Leader(Vector3 targetP)
    {
        //have the actors follow a leader as a group
        
        //if this actor is the leader, have the actor seek the target given
        if (rigidbody == targetScript.agents[0].rigidbody)
            return Seek(targetP);
            
        //otherwise, have the actor arrive at the leader's position while keeping some
        //distance from its neighbors
        Vector3 leaderFollow = 5*Separation()+5*Arrival(targetScript.agents[0].rigidbody.position);
        
        //Have the actor move at max speed in the calculated direction
        return g_fMaxSpeed*leaderFollow.normalized;
    }
    
    Vector3 Seek(Vector3 targetP)
    {
        //We want the actor to go as fast as he can toward the target
        //Vd = maxspeed * normalized_e
        
        //e = targetPosition - currentPosition
        Vector3 e = targetP-rigidbody.position;
        
        //make sure actor is not already at target
        float dist = e.magnitude;
        //if he is, return the zero vector.  actor need not move.
        if (dist == 0)
            return Vector3.zero;
        
        //e = ||e|| * maxspeed
        e = e.normalized * g_fMaxSpeed;
        
        return e;   
    }
    
    Vector3 Separation()
    {  
        //Have the actors keep some distance from each other.
        //Each actor tries to separate a little from its neighbors.
        
        Vector3 sumDist = Vector3.zero;
        Vector3 dist = Vector3.zero;
        
        //Add up the distance over distance length squared for every neighbor within a radius
        for (int i = 0; i <  targetScript.agentNumber ; i++)
        {
            dist = rigidbody.position-targetScript.agents[i].rigidbody.position;
            if (dist.magnitude < g_fKNeighborhood && rigidbody != targetScript.agents[i].rigidbody)
            {
                sumDist+=dist/(dist.magnitude*dist.magnitude);
            }
        }
        
        //return this sum weighted by g_fKSeparation
        return g_fKSeparation*(sumDist);
    }
    
    
    Vector3 Wander()
    {
        //Have the actor wander the board aimlessly.
        
        //Come up with a noise vector based on the actor's current velocity.
        Vector3 rnoise = g_fKNoise*(new Vector3(Random.Range(-rigidbody.velocity.magnitude, rigidbody.velocity.magnitude), 0, Random.Range(-rigidbody.velocity.magnitude, rigidbody.velocity.magnitude)).normalized);
        
        //Add this noise to the actor's current wandering velocity addition factor to get
        //a new wandering velocity addition factor.
        m_vWander = g_fKWander*((m_vWander+rnoise));
        
        //Add the wandering velocity addition factor to the current velocity.
        //Have the actor move at max speed in the calculated direction.
        return g_fMaxSpeed*((rigidbody.velocity+m_vWander).normalized);
    }
}