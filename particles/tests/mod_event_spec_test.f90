!> This module contains testcases for the event system
module mod_event_spec_test
use fruit
use mod_event
use mod_particle_sim
implicit none
private
public :: run_fruit_event_spec
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_event_spec
  implicit none
  write(*,'(/A)') "  ... setting-up: event spec"
  write(*,'(/A)') "  ... running: event spec"
  call run_test_case(test_create_event_with_stop_action,'test_create_event_with_stop_action')
  call run_test_case(test_event_run_at,'test_event_run_at')
  call run_test_case(test_next_event_at,'test_next_event_at')
  write(*,'(/A)') "  ... tearing-down: event spec"
end subroutine run_fruit_event_spec

!> Tests ------------------------------------------

subroutine test_create_event_with_stop_action
  type(event), dimension(:), allocatable :: events
  events = [event(stop_action(), start=1.d0)]
  call assert_equals(events(1)%start, 1.d0)
end subroutine test_create_event_with_stop_action

subroutine test_event_run_at
  type(event) :: e
  e = event(stop_action(), start=1.d0)
  call assert_false(e%run_at(0.d0), 'must not run before start')
  call assert_true(e%run_at(1.d0), 'must run at start')
  call assert_false(e%run_at(2.d0), 'must not run after start if no step is set')
  e = event(stop_action(), start=1.d0, step=2.d0)
  call assert_true(e%run_at(3.d0), 'must run at start + 1 step')
  call assert_false(e%run_at(3.1d0), 'must not run at start + 1.05 step')
  call assert_false(e%run_at(3.0001d0), 'must not run at start + 1.00005 step')
  call assert_false(e%run_at(-1.d0), 'must not run at start - 1 step')
  e = event(stop_action(), start=1.d0, step=2.d0, end=3.d0)
  call assert_true(e%run_at(3.d0), 'must run at start + 1 step == end')
  call assert_true(e%run_at(3.d0+tick/2), 'must run at start + 1 step == end')
  call assert_false(e%run_at(5.d0), 'must not run > end')
  e = event(stop_action(), start=1.d0, end=2.d0)
  call assert_true(e%run_at(1.d0), 'must run at start')
  call assert_false(e%run_at(2.d0), 'must not run after start if no step is set')
end subroutine test_event_run_at

subroutine test_next_event_at
  type(particle_sim) :: sim
  type(event), dimension(:), allocatable :: events
  events = [event(stop_action(), start=1.d0, end=4.d0, step=0.1d0)]
  sim%time = 0.d0;   call assert_equals(1.d0, next_event_at(sim, events), tick, 'must find start if in the future')
  sim%time = 1.05d0; call assert_equals(1.1d0, next_event_at(sim, events), tick, 'must find with step if in the future')
  sim%time = 1.1d0;  call assert_equals(1.2d0, next_event_at(sim, events), tick, 'must skip events at this exact time')
  sim%time = 3.95d0
  call assert_equals(4.0d0, next_event_at(sim, events), tick, 'must return the time at the last event')
  call assert_true(sim%stop_now, 'must set stop_now just at the last event')
  sim%time = 4.d0
  sim%stop_now = .false.
  call assert_equals(4.0d0, next_event_at(sim, events), tick, 'must return the time at the last event')
  call assert_true(sim%stop_now, 'must set stop_now at the last event')
  sim%time = 5.d0
  sim%stop_now = .false.
  call assert_equals(5.0d0, next_event_at(sim, events), tick, 'must return the current time if past the last event')
  call assert_true(sim%stop_now, 'must set stop_now past the last event')

  events = [event(stop_action(), start=1.d0, end=4.d0, step=0.1d0), &
            event(stop_action(), start=5.d0)]
  sim%time = 3.95d0; call assert_equals(4.d0, next_event_at(sim, events), tick, 'must find step if in region')
  sim%time = 4.d0;   call assert_equals(5.d0, next_event_at(sim, events), tick, 'must find start even after end of other')
  sim%time = 5.d0;   call assert_equals(5.d0, next_event_at(sim, events), tick, 'must find start even after end of other')
  call assert_true(sim%stop_now, 'must set stop_now just at the last event')
  sim%time = 6.d0; sim%stop_now = .false.
  call assert_equals(6.d0, next_event_at(sim, events), tick, 'must find start even after end of other')
  call assert_true(sim%stop_now, 'must set stop_now after the last event')
end subroutine test_next_event_at

!> ------------------------------------------------
end module mod_event_spec_test
