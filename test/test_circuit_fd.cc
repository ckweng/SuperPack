#include <catch2/catch.hpp>
#include <iostream>

#include "tp/circuits.h"

#define PARTY for(std::size_t i = 0; i < n_parties; i++)
#define BATCH for(std::size_t i = 0; i < batch_size; i++)

TEST_CASE("Secure computation function-dependent") {
  SECTION("Hand-crafted circuit") {
    // Inputs x,y,u,v
    // x' = (x+y)*x,  y' = (x+y)*y
    // u' = (u+v)*u,  v' = (u+v)*v
    // z = (x' + y') + (u' + v')

    std::size_t threshold = 4; // has to be even
    std::size_t batch_size = (threshold + 2)/2;
    std::size_t n_parties = threshold + 2*(batch_size - 1) + 1;
    auto networks = scl::Network::CreateFullInMemory(n_parties);
    std::size_t n_clients = n_parties;

    tp::FF lambda(649823289);

    std::vector<tp::Circuit> circuits;
    circuits.reserve(n_parties);

    PARTY {
      auto c = tp::Circuit(n_clients, batch_size);
      
      // input batches x, y, u, v
      auto x = c.Input(0);
      auto y = c.Input(1);
      auto u = c.Input(0);
      auto v = c.Input(1);

      c.CloseInputs();
      
      std::vector<std::shared_ptr<tp::Gate>> xPy;
      std::vector<std::shared_ptr<tp::Gate>> uPv;
      BATCH {
        xPy.push_back(c.Add(x->GetInputGate(i), y->GetInputGate(i)));
        uPv.push_back(c.Add(u->GetInputGate(i), v->GetInputGate(i)));
      }

      std::vector<std::shared_ptr<tp::Gate>> x_;
      std::vector<std::shared_ptr<tp::Gate>> y_;
      std::vector<std::shared_ptr<tp::Gate>> u_;
      std::vector<std::shared_ptr<tp::Gate>> v_;
      BATCH {
        x_.push_back(c.Mult(xPy[i], x->GetInputGate(i)));
        y_.push_back(c.Mult(xPy[i], y->GetInputGate(i)));
        u_.push_back(c.Mult(uPv[i], u->GetInputGate(i)));
        v_.push_back(c.Mult(uPv[i], v->GetInputGate(i)));
      }
      c.LastLayer();
   
      std::vector<std::shared_ptr<tp::Gate>> z1;
      std::vector<std::shared_ptr<tp::Gate>> z2;
      BATCH {
        z1.push_back(c.Add(x_[i], y_[i]));
        z2.push_back(c.Add(u_[i], v_[i]));
      }

      std::vector<std::shared_ptr<tp::Gate>> z;
      BATCH {
        z.push_back(c.Add(z1[i], z2[i]));
      }

      auto output = c.Output(0,z);
      c.CloseOutputs();

      c.SetNetwork(std::make_shared<scl::Network>(networks[i]), i);
      c.GenCorrelator();
      c._DummyPrepFI(lambda);
      c.MapCorrToCircuit();

      circuits.emplace_back(c);
    }

    // circuit-dependent preprocessing
    PARTY { circuits[i].PrepMultPartiesSendP1(); }
    PARTY { circuits[i].PrepMultP1Receives(); }
    PARTY { circuits[i].PrepIO(); }

    // SET INPUTS
    std::vector<tp::FF> X(batch_size, tp::FF(21321));
    std::vector<tp::FF> Y(batch_size, tp::FF(-3421));
    std::vector<tp::FF> U(batch_size, tp::FF(170942));
    std::vector<tp::FF> V(batch_size, tp::FF(-894));

    // P1 sets inpus X and U
    circuits[0].SetInputs(std::vector<std::vector<tp::FF>>{X, U});
    // P2 sets inpus Y and V
    circuits[1].SetInputs(std::vector<std::vector<tp::FF>>{Y, V});

    // Input protocol
    REQUIRE(circuits[0].GetInputGate(0,0,0)->IsLearned() == false);
    REQUIRE(circuits[0].GetInputGate(0,1,0)->IsLearned() == false);
    REQUIRE(circuits[0].GetInputGate(1,0,0)->IsLearned() == false);
    REQUIRE(circuits[0].GetInputGate(1,1,0)->IsLearned() == false);

    PARTY { circuits[i].InputPartiesRevealShares(); }
    PARTY { circuits[i].InputClientReconstructAndDistributes(); }
    PARTY { circuits[i].InputPartiesReceive(); }
    PARTY { circuits[i].InputClientSendMuToP1(); }
    PARTY { circuits[i].InputP1ReceiveMu(); }
    PARTY { circuits[i].InputPartiesAuthenticateMu(); }

    BATCH {
      REQUIRE(circuits[0].GetInputGate(0,0,i)->GetMu() == X[i] - lambda);
      REQUIRE(circuits[0].GetInputGate(0,1,i)->GetMu() == U[i] - lambda);
      REQUIRE(circuits[0].GetInputGate(1,0,i)->GetMu() == Y[i] - lambda);
      REQUIRE(circuits[0].GetInputGate(1,1,i)->GetMu() == V[i] - lambda);
    }

    // Multiplications (there is only one layer)
    for(size_t i = 0; i < 4*batch_size; ++i) {
      REQUIRE(circuits[0].GetMultGate(0,i)->IsLearned() == false);
    }
    
    PARTY { circuits[i].MultP1Distribute(0); }
    PARTY { circuits[i].MultPartiesReceive(0); }
    PARTY { circuits[i].MultPartiesMultiplyAndSend(0); }
    PARTY { circuits[i].MultP1ReceivesAndReconstruct(0); }
    PARTY { circuits[i].MultPartiesCheckAndComputeAuthMu(0); }

    for(size_t i = 0; i < batch_size; ++i) {
      REQUIRE(circuits[0].GetMultGate(0,4*i)->GetMu() == (X[i]+Y[i])*X[i] - lambda);
      REQUIRE(circuits[0].GetMultGate(0,4*i+1)->GetMu() == (X[i]+Y[i])*Y[i] - lambda);
      REQUIRE(circuits[0].GetMultGate(0,4*i+2)->GetMu() == (U[i]+V[i])*U[i] - lambda);
      REQUIRE(circuits[0].GetMultGate(0,4*i+3)->GetMu() == (U[i]+V[i])*V[i] - lambda);
    }
    
    // Output protocol
    PARTY { circuits[i].OutputPartiesRevealShares(); }
    PARTY { circuits[i].OutputP1ReconstructAndDistributes(); }
    PARTY { circuits[i].OutputPartiesReceiveShares(); }
    PARTY { circuits[i].OutputPartiesVerify(); }
    PARTY { circuits[i].OutputPartiesCheckTheta(); }
    PARTY { circuits[i].OutputPartiesRevealOutputVal(); }
    PARTY { circuits[i].OutputClientReconstructSecret(); }


    // Check output
    tp::FF real = (X[0]+Y[0])*X[0] + (X[0]+Y[0])*Y[0] + (U[0]+V[0])*U[0] + (U[0]+V[0])*V[0];
    BATCH { REQUIRE(circuits[0].GetOutputGate(0,0,i)->GetValue() == real); }
  }
  
  /*SECTION("Artificial circuit") 
    {
      std::size_t batch_size = 2;
      std::size_t n_parties = 4*batch_size - 3;

      tp::CircuitConfig config;
      config.n_parties = n_parties;
      config.inp_gates = std::vector<std::size_t>(n_parties, 0);
      config.inp_gates[0] = 2;
      config.out_gates = std::vector<std::size_t>(n_parties, 0);
      config.out_gates[0] = 2;
      config.width = 100;
      config.depth = 2;
      config.batch_size = batch_size;
	  
      auto networks = scl::Network::CreateFullInMemory(n_parties);

      std::vector<tp::Circuit> circuits;
      circuits.reserve(n_parties);

      tp::FF lambda(-45298432);

      PARTY {
	auto c = tp::Circuit::FromConfig(config);
	c.SetNetwork(std::make_shared<scl::Network>(networks[i]), i);
	c._DummyPrep(lambda);
	circuits.emplace_back(c);
      }

      std::vector<tp::FF> inputs{tp::FF(0432432), tp::FF(54982)};
      circuits[0].SetClearInputsFlat(inputs);
      auto result = circuits[0].GetClearOutputsFlat();
      circuits[0].SetInputs(inputs);

      // INPUT
      PARTY { circuits[i].InputOwnerSendsP1(); }
      PARTY { circuits[i].InputP1Receives(); }

      // MULT
      for (std::size_t layer = 0; layer < config.depth; layer++) {
	PARTY { circuits[i].MultP1Sends(layer); }
	PARTY { circuits[i].MultPartiesReceive(layer); }
	PARTY { circuits[i].MultPartiesSend(layer); }
	PARTY { circuits[i].MultP1Receives(layer); }
      }

      // OUTPUT
      PARTY { circuits[i].OutputP1SendsMu(); }
      PARTY { circuits[i].OutputOwnerReceivesMu(); }

      // Check output
      REQUIRE(circuits[0].GetOutputs() == result);
    }*/
}
