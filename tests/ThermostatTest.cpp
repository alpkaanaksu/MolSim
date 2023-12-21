//
// Created by Feryal Ezgi on 21.12.2023.
//
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

#include "../src/models/LinkedCellParticleContainer.h"
#include "../src/Simulation.h"

class ThermostatTest : public ::testing::Test {
protected:
    bool notNear(double actual, double expected, double epsilon) {
        return std::abs(actual - expected) > epsilon;
    }

};

// Verifies that the temperature is not adjusted if a thermostat is not given
TEST_F(ThermostatTest, NoThermostat) {
    std::string in = "../input/test/thermostat/no_thermostat.json";
    Simulation simulation(in);

    auto particles = simulation.getParticles();
    auto linkedCellParticles = std::dynamic_pointer_cast<LinkedCellParticleContainer>(particles);
    ASSERT_NE(linkedCellParticles, nullptr);

    simulation.run();
    EXPECT_NEAR(simulation.getThermostat().getInitialTemperature(), 0, 1e-3);
    EXPECT_NEAR(simulation.getThermostat().getTargetTemperature(), 0, 1e-3);

    // The temperature after the simulation shouldn't be near 20 degrees as per the input molecule velocities, gravity etc.
    EXPECT_TRUE(notNear(simulation.getThermostat().getCurrentTemperature(*particles), 20, 1e-3));
}


TEST_F(ThermostatTest, Heating15To20) {
    std::string in = "../input/test/thermostat/heating_15_20.json";
    Simulation simulation(in);

    auto particles = simulation.getParticles();
    auto linkedCellParticles = std::dynamic_pointer_cast<LinkedCellParticleContainer>(particles);
    ASSERT_NE(linkedCellParticles, nullptr);

    simulation.run();
    SPDLOG_INFO("Current Temperature: {}", simulation.getThermostat().getCurrentTemperature(*particles));

    EXPECT_EQ(simulation.getThermostat().getInitialTemperature(), 15);
    EXPECT_EQ(simulation.getThermostat().getTargetTemperature(), 20);
    EXPECT_NEAR(simulation.getThermostat().getCurrentTemperature(*particles), 20, 1e-3);

}

TEST_F(ThermostatTest, Heating0To25) {
    std::string in = "../input/test/thermostat/heating_0_25.json";
    Simulation simulation(in);

    auto particles = simulation.getParticles();
    auto linkedCellParticles = std::dynamic_pointer_cast<LinkedCellParticleContainer>(particles);
    ASSERT_NE(linkedCellParticles, nullptr);

    simulation.run();
    SPDLOG_INFO("Current Temperature: {}", simulation.getThermostat().getCurrentTemperature(*particles));

    EXPECT_EQ(simulation.getThermostat().getInitialTemperature(), 0);
    EXPECT_EQ(simulation.getThermostat().getTargetTemperature(), 25);
    EXPECT_NEAR(simulation.getThermostat().getCurrentTemperature(*particles), 25, 1e-3);

}

TEST_F(ThermostatTest, Heating3D0To25) {
    std::string in = "../input/test/thermostat/heating_3D_0_25.json";
    Simulation simulation(in);

    auto particles = simulation.getParticles();
    auto linkedCellParticles = std::dynamic_pointer_cast<LinkedCellParticleContainer>(particles);
    ASSERT_NE(linkedCellParticles, nullptr);

    simulation.run();
    SPDLOG_INFO("Current Temperature: {}", simulation.getThermostat().getCurrentTemperature(*particles));

    EXPECT_EQ(simulation.getThermostat().getInitialTemperature(), 0);
    EXPECT_EQ(simulation.getThermostat().getTargetTemperature(), 25);
    EXPECT_NEAR(simulation.getThermostat().getCurrentTemperature(*particles), 25, 1e-3);

}

TEST_F(ThermostatTest, HeatingGradual2To6) {
    std::string in = "../input/test/thermostat/heating_gradual_2_6.json";
    Simulation simulation(in);

    auto particles = simulation.getParticles();
    auto linkedCellParticles = std::dynamic_pointer_cast<LinkedCellParticleContainer>(particles);
    ASSERT_NE(linkedCellParticles, nullptr);

    simulation.run();
    SPDLOG_INFO("Current Temperature: {}", simulation.getThermostat().getCurrentTemperature(*particles));

    EXPECT_EQ(simulation.getThermostat().getInitialTemperature(), 2);
    EXPECT_EQ(simulation.getThermostat().getTargetTemperature(), 6);
    EXPECT_NEAR(simulation.getThermostat().getCurrentTemperature(*particles), 6, 1e-3);
}

TEST_F(ThermostatTest, HeatingGradualSmallMaxTemperature) {
    std::string in = "../input/test/thermostat/heating_gradual_smallMaxTemp.json";
    Simulation simulation(in);

    auto particles = simulation.getParticles();
    auto linkedCellParticles = std::dynamic_pointer_cast<LinkedCellParticleContainer>(particles);
    ASSERT_NE(linkedCellParticles, nullptr);

    simulation.run();
    SPDLOG_INFO("Current Temperature: {}", simulation.getThermostat().getCurrentTemperature(*particles));

    EXPECT_EQ(simulation.getThermostat().getInitialTemperature(), 2);
    EXPECT_EQ(simulation.getThermostat().getTargetTemperature(), 20);

    EXPECT_GE(simulation.getThermostat().getCurrentTemperature(*particles),
              simulation.getThermostat().getInitialTemperature());
    EXPECT_LE(simulation.getThermostat().getCurrentTemperature(*particles),
              simulation.getThermostat().getTargetTemperature());

    EXPECT_TRUE(notNear(simulation.getThermostat().getCurrentTemperature(*particles), 20, 1e-3));
}
